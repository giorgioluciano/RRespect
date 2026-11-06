# ReSpectR

An R package for continuous and discrete relaxation spectrum extraction from time-domain relaxation modulus data and frequency-domain storage and loss moduli. Python is not required to run the R package.

## Install the validated reference

```r
# install.packages("remotes")
remotes::install_github("giorgioluciano/RRespect@v0.5.0")
library(ReSpectR)
```

The historical validation reference is tag `v0.5.0`, commit `16291f92f2f8c0034c554af9866e65ea555c4f0c`. The `submission-cleanup` branch is a work in progress, not a new validated release.

## Basic workflow

```r
library(ReSpectR)
par <- setParams(
  domain = "time", dataFile = "Gt.dat",
  solver = "experimental",
  verbose = FALSE, plotting = FALSE
)
crs <- getContinuousSpectrum(par, writeOutput = FALSE)
drs <- getDiscreteSpectrum(par, crs = crs, writeOutput = FALSE)
```

For frequency-domain data, use `domain = "frequency"` and a file containing angular frequency, storage modulus and loss modulus. See `?setParams` for input formats and options.

The updated default workflow uses trust-region reflective nonlinear least squares. The selector `solver = "experimental"` is retained for this workflow; historical solvers remain available through `solver = "legacy"`. Continuous spectra use regularisation with Bayesian parameter selection. Discrete spectra use non-negative fitting, refinement, pruning and mode merging.

## Documentation and validation

- Package help: `help(package = "ReSpectR")`.
- Vignettes, when installed: `browseVignettes("ReSpectR")`.
- Validation record: [inst/VALIDATION-2026-10-07.md](inst/VALIDATION-2026-10-07.md).
- General validation notes: [inst/VALIDATION.md](inst/VALIDATION.md).

The recorded Python comparison covered 23 configurations: 20 passed all applicable recorded numerical criteria and 22 returned matching final mode counts. Discrepancies remain in frequency test1u, frequency test3 and time test2. These results refer to the documented reference revisions and settings, not arbitrary subsequent changes. Agreement between implementations does not establish uniqueness or physical correctness of a recovered spectrum, or a good fit to every dataset.

## Development

From a local checkout with development dependencies installed:

```r
devtools::test()
devtools::check(cran = TRUE, manual = FALSE)
```

Historical prototypes and generated artifacts under `legacy/` have been removed from this branch. They remain available in the historical `v0.5.0` tag. Required package data remain in `inst/extdata`, and the `solver = "legacy"` implementation under `R/` is retained. Upstream notices are preserved in `inst/COPYRIGHTS`. A separate snapshot of commit `fe79ff849ae62b8effe35385545a55d646b7b16e`, with `legacy/` removed, completed R CMD check with 0 errors, 0 warnings and 2 notes (new submission; current time could not be verified). The launching session used R 4.5.2 on Windows 11. This is not a claim that the subsequent final cleanup commit has already passed new checks; final package checks remain required.

## Licence and support

The package currently declares GPL-3. Upstream notices and attribution must be reviewed before final submission; this cleanup does not relicense third-party material.

Report issues at https://github.com/giorgioluciano/RRespect/issues.
