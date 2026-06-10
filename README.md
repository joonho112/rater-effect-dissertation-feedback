# Methodological Notes on Two Rater-Effect Simulation Studies

Methodological feedback on a doctoral student's dissertation — research notes on the two
rater-effect simulation studies at its methodological core (a Many-Facet Rasch Model study and
a Hierarchical Rater Model study), worked toward a journal manuscript. Built as a
[Quarto](https://quarto.org) book.

**Author:** JoonHo Lee (jlee296@ua.edu) · University of Alabama

## About

The notes take each methodological concern through a four-part path — what the dissertation
does, why it threatens the conclusion, what would fix it, and what a reviewer will ask —
across seven chapters, an executive summary, and a revision roadmap. Where a numerical claim
is checked, the relevant equations are re-implemented from scratch in R (see `R/common.R`),
independent of the original estimation software; the figures are reproducible by rendering the
book.

## Read online

The rendered book is published with GitHub Pages:
<https://joonho112.github.io/rater-effect-dissertation-feedback/>

## Build locally

```bash
quarto render     # render the full book to _book/
quarto preview    # live preview while editing
```

Requires [Quarto](https://quarto.org) and R with the packages `ggplot2`, `patchwork`, and
`statmod`.

## Repository layout

```
index.qmd                Preface
00-…-08-….qmd            Executive summary, seven chapters, revision roadmap
99-reproducibility.qmd   Reproducibility appendix
R/common.R               Shared equations (MFRM / GPCM / signal-detection), in R
_quarto.yml              Book configuration
```

## License

Released under the MIT License (see [`LICENSE`](LICENSE)).
