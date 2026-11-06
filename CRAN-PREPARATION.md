# GPcluster release preparation

Prepared November 4, 2026. Changes remain local; not submitted to CRAN.

## Implemented

- Native mlegp nugget estimation; laGP dependency removed.
- Selection among candidate K values by minimum existing LOOCV RMSE.
- K = 1 baseline and size-constrained initialization (minimum 3*d).
- Log-density E-step assignments, freshly computed weights, and correct handling
  of staying in a cluster already at its maximum size.
- Schur-complement inverse updates, with full Cholesky inversion for nearly
  singular additions. Updated model matrices are returned from the E-step.
- Adaptive mixture quantile bounds and zero-variance component handling.
- Input validation, runnable examples, regression tests, and metadata cleanup.
- Published Statistica Sinica citation (2023), 33(2), 893–918,
  doi:10.5705/ss.202020.0456; available through citation("GPcluster").
- GitHub Actions for Windows/macOS/Linux R-release and Linux R-devel, plus
  a separate two-worker PSOCK integration test.

## Evidence

Local R 4.4.3 / Apple Silicon macOS: source-package check including PDF manual
passed with 0 errors and 0 warnings, with incoming checks disabled. Notes were
unverifiable current time and an older HTML validator rejecting HTML5 generated
by R. Logs are under /tmp/gpcluster-release-check/GPcluster.Rcheck.

Two-worker PSOCK fitting/prediction passed locally for K = 1 and K = 2 with
estimated nuggets. Regression tests compare covariance updates with full inverses,
including near-duplicate inputs and actual E-step transfers. They also cover
selection, nugget estimation, quantiles, invalid inputs, and initialization.

The LOOCV score conditions on fitted partitions/parameters, rather than refitting
all clustering for every omitted point. Use held-out data for final performance.

## Still required for release

Run the prepared workflow after pushing the changes to GitHub and inspect actual
Windows/current-release/devel results. Local R 4.4.3 does not establish those.
Review online incoming/URL results, confirm maintainer address and authorship,
and prepare final submission comments. No CRAN submission has been made.

Policy: https://cran.r-project.org/web/packages/policies.html

Online check completed with 0 errors, 0 warnings, and 4 notes: new submission,
current-time verification, missing pandoc in that environment, and the legacy HTML
validator. The subsequent one-dimensional initialization shape fix passed the
final local check (0 errors, 0 warnings, 2 notes) and all regression tests.
Saved logs: ../GPcluster-local-check.log and ../GPcluster-online-check.log.
Final source archive: ../GPcluster_0.1.tar.gz.
