# ReSpectR 0.5.0: validation record

Date: 2026-10-07.

## Package checks

- Updated TRF workflow is the default (`solver = "experimental"`).
- Historical workflow remains available with `solver = "legacy"`.
- Full testthat suite completed without reported failures.
- Local R CMD check completed with Status: OK.
- Check options: --no-manual --as-cran.
- Examples, tests, vignette building and rebuilding passed.
- Environment: R 4.6.0, Windows 11, x86_64-w64-mingw32.
- This record does not imply CRAN submission or acceptance.

## Python implementation comparison

References:
- Time: cc2545c461deda614d0344c0731d9d5e55c04f43.
- Frequency: 908559720322941622ab39ced82e376776a39183.
- NumPy 2.5.3; SciPy 1.18.1.

The comparison covered 23 configurations:
8 upstream time-domain configurations, 11 upstream frequency-domain
configurations, and 4 synthetic forced-merge configurations.

20 configurations passed all recorded comparison criteria.
3 configurations exhibited discrepancies:

- frequency-test1u: final mode counts agreed (17), but coefficients,
  times, uncertainties and fitted curves exceeded comparison thresholds.
  Maximum curve difference normalised by observations: 1.53e-4.
- frequency-test3: final mode counts differed (R: 14; Python: 12).
  Weighted SSE: approximately 0.03342427 and 0.02682667.
  Maximum observation-normalised curve difference: approximately 3%.
- time-test2: final mode counts agreed (4); coefficients, times and
  fitted curves passed, but merge counts and uncertainties differed.
  Maximum observation-normalised curve difference: 6.76e-8.

All four forced-merge configurations passed and each executed one merge.
All recorded continuous-spectrum and AIC-scan controls passed
their stated comparison thresholds.

## Numerical sensitivity

Additional frequency-test3 diagnostics showed that both implementations
can reject a refined solution after exp(log(tau)) moves a feasible
boundary mode just outside the permitted window.
Python also rejected refinement from the original R inputs exported
with verified double-precision round-trip preservation.
Internal optimisation paths and pruned mode counts were not identical.

This is evidence of shared boundary sensitivity, not proof of complete
end-to-end equivalence or physical superiority of either implementation.
The original three comparison discrepancies remain documented.

## Reproducibility settings

Spectrum grid: ns = 100; freqEnd = 1; smFacLam = 0.
Automatic lambda range: 1e-10 to 1e3; density = 2.
Baseline-weight step: 0.2; no imposed mode-count cap.
Ordinary minimum tau spacing: 1.25.
Forced-merge cases: lambda = 0.01; minimum tau spacing = 30.
Upstream configurations used 100 observation points;
weighted synthetic configurations retained 30 observation points.

Counts and merge counts were compared exactly.
Curve tolerances: rtol = 1e-5, atol = 1e-9.
Discrete coefficients/times: rtol = 1e-3, atol = 1e-9.
Discrete time uncertainties: rtol = 1e-2, atol = 1e-8.
The detailed comparison scripts and CSV reports contain the remaining
per-field thresholds and diagnostics.

The detailed audit artifacts are retained locally under output/;
they are not automatically included in the package or Git commit.
