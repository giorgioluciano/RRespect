# Pre-submission validation

## Package quality
- [ ] Clean installation from source
- [x] `devtools::test()` passes — 40 passed, 0 failed, 0 test warnings on Windows 11 / R 4.6.0 (2026-10-05)
- [x] Test files run individually with no global-environment pollution
- [ ] `devtools::check()` passes
- [x] `R CMD check --as-cran` passes — 0 errors, 0 warnings; CRAN incoming notes reviewed on Windows 11 / R 4.6.0 (2026-10-05)
- [x] PDF and HTML package manuals build successfully
- [x] Package vignettes build and rebuild successfully
- [ ] Examples run
- [ ] Vignettes build

## Numerical validation
- [ ] Synthetic one-mode Maxwell case
- [ ] Synthetic multi-mode Maxwell case
- [ ] Noisy synthetic case
- [ ] Boundary and invalid-input cases
- [ ] Regression fixtures frozen
- [ ] Results reproduced from manuscript scripts

## Release
- [ ] LICENSE selected and committed
- [ ] README installation instructions tested
- [ ] CITATION file present
- [ ] NEWS.md updated
- [ ] Version bumped
- [ ] GitHub Actions green
- [ ] Release tag created
- [ ] Independent co-author verification complete
