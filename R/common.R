# =============================================================================
# R/common.R — shared helpers for "Methodological Notes on Two Rater-Effect
# Simulation Studies".
#
# These are a faithful R port of an independent Python reference implementation. They
# implement the dissertation's OWN equations from scratch (they do NOT call TAM
# or immer), so that each critique can be reproduced numerically:
#
#   * PC-MFRM data-generating model   -> Eq. (7)/(8)   pcmfrm_probs()
#   * MFRM rater fit (infit/outfit)   -> Eq. (2)/(3)   rater_fit()
#   * GPCM (HRM stage 1)              -> Eq. (5)        gpcm_probs()
#   * Signal Detection Model (stage 2)-> Eq. (6)        sdm_probs()
#
# Scale: 5 ordered categories scored 0..4 (K = 5), four fixed thresholds
#   tau = (-1.5, -0.5, 0.5, 1.5).
#
# This file was cross-checked against an independent Python implementation of the same
# equations; deterministic functions agree to < 1e-12.
# =============================================================================

# ---- Constants --------------------------------------------------------------
K    <- 5L
TAU  <- c(-1.5, -0.5, 0.5, 1.5)   # dissertation's fixed thresholds
CATS <- 0:4                        # category scores 0,1,2,3,4

# ---- Colour palette (mirrors the Python C_* constants) ----------------------
pal_companion <- c(
  none   = "#4C4C4C",
  severe = "#C0392B",
  misfit = "#2471A3",
  accent = "#117A65",
  warn   = "#B9770E"
)

# ---- internal: row-wise log-sum-exp stabilisation ---------------------------
# Subtract each row's maximum (softmax is shift-invariant per row, so this only
# improves numerical stability and never changes the normalised result).
.lse_stabilise <- function(m) sweep(m, 1L, apply(m, 1L, max), "-")

# ---- PC-MFRM (Eq. 7/8) ------------------------------------------------------
#' Category probabilities P(X=k | theta) for one (rater, item).
#'   ln[ P(x=k)/P(x=k-1) ] = alpha * (theta - lam - delta - tau_k)
#' alpha = 1 -> standard PC-MFRM ("accurate" rater);
#' alpha = 0 -> every adjacent log-odds is 0 -> UNIFORM over categories.
#' @return N x K matrix (rows = persons, cols = categories 0..4); 1xK if scalar.
pcmfrm_probs <- function(theta, lam = 0, delta = 0, tau = TAU, alpha = 1) {
  theta   <- as.numeric(theta)
  N       <- length(theta)
  tau_cum <- c(0, cumsum(tau))                 # length K
  base    <- theta - lam - delta               # length N
  # cum[i,k] = alpha * (CATS[k]*base[i] - tau_cum[k])
  cum <- alpha * (outer(base, as.numeric(CATS)) - outer(rep(1, N), tau_cum))
  cum <- .lse_stabilise(cum)
  w   <- exp(cum)
  matrix(w / rowSums(w), nrow = N, ncol = K)
}

#' Expected score E[X] from a probability matrix.
expected_score <- function(probs) as.numeric(probs %*% as.numeric(CATS))

#' Score variance Var[X] = E[X^2] - E[X]^2 (the model information weight W).
score_variance <- function(probs) {
  ex  <- expected_score(probs)
  ex2 <- as.numeric(probs %*% (as.numeric(CATS)^2))
  ex2 - ex^2
}

#' Draw one integer rating 0..K-1 per row via inverse-CDF on runif().
#' Caller sets the RNG seed. Mirrors Python (u[:,None] > cdf).sum(axis=1).
simulate_rating <- function(probs) {
  if (is.null(dim(probs))) probs <- matrix(probs, nrow = 1)
  N   <- nrow(probs)
  u   <- runif(N)
  cdf <- t(apply(probs, 1L, cumsum))           # N x K
  as.integer(rowSums(matrix(u, nrow = N, ncol = K) > cdf))
}

#' MFRM rater fit statistics (Eq. 2 outfit, Eq. 3 infit).
#' Standardised residual Z = (X - E[X]) / sqrt(Var[X]); W = Var[X].
#' @return list(infit, outfit). Variance clipped at 1e-9.
rater_fit <- function(observed, probs) {
  ex  <- expected_score(probs)
  vr  <- pmax(score_variance(probs), 1e-9)
  z2  <- (as.numeric(observed) - ex)^2 / vr
  list(infit = sum(z2 * vr) / sum(vr), outfit = mean(z2))
}

# ---- GPCM (HRM stage 1, Eq. 5) ----------------------------------------------
#' GPCM ideal-rating probabilities.
#'   P(xi=k|theta) propto exp( sum_{j<=k} a*(theta - b) - gamma_j )
#' Default gamma = cumsum(TAU).
gpcm_probs <- function(theta, a = 1, b = 0, gamma = NULL) {
  if (is.null(gamma)) gamma <- cumsum(TAU)
  theta   <- as.numeric(theta)
  N       <- length(theta)
  gam_cum <- c(0, gamma)                        # length K
  base    <- a * (theta - b)                    # length N
  cum <- outer(base, as.numeric(CATS)) - outer(rep(1, N), gam_cum)
  cum <- .lse_stabilise(cum)
  w   <- exp(cum)
  matrix(w / rowSums(w), nrow = N, ncol = K)
}

# ---- Signal Detection Model (HRM stage 2, Eq. 6) ----------------------------
#' SDM category probabilities for ideal rating(s) xi.
#'   p(X=k | xi) propto exp( -1/(2 psi^2) * (k - (xi + phi))^2 )
#' As psi -> Inf the distribution -> UNIFORM(1/K) for ANY xi.
sdm_probs <- function(xi, phi = 0, psi = 0.4) {
  xi     <- as.numeric(xi)
  N      <- length(xi)
  centre <- xi + phi                            # length N
  d2  <- (outer(rep(1, N), as.numeric(CATS)) - outer(centre, rep(1, K)))^2
  cum <- -d2 / (2 * psi^2)
  cum <- .lse_stabilise(cum)
  w   <- exp(cum)
  matrix(w / rowSums(w), nrow = N, ncol = K)
}

#' Total log-likelihood of a rater's ratings under the SDM.
#' observed_X are integers 0..K-1; clamped, then column = value + 1.
sdm_loglik <- function(observed_X, xi, phi, psi) {
  p   <- sdm_probs(xi, phi = phi, psi = psi)    # M x K
  M   <- nrow(p)
  idx <- pmax(pmin(as.integer(observed_X), K - 1L), 0L)
  pk  <- p[cbind(seq_len(M), idx + 1L)]
  sum(log(pmax(pk, 1e-300)))
}

# ---- ggplot2 theme (loaded lazily so sourcing never needs ggplot2) ----------
#' Minimal, clean theme for companion figures. Requires ggplot2 at call time.
theme_companion <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for theme_companion(); please install it.")
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major    = ggplot2::element_line(colour = "grey92", linewidth = 0.4),
      panel.grid.minor    = ggplot2::element_blank(),
      axis.line.x         = ggplot2::element_line(colour = "grey60", linewidth = 0.4),
      axis.line.y         = ggplot2::element_line(colour = "grey60", linewidth = 0.4),
      legend.position     = "bottom",
      plot.title          = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption        = ggplot2::element_text(hjust = 0, colour = "grey30")
    )
}
