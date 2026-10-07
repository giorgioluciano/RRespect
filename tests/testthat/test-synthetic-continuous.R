test_that("a smooth continuous spectrum is recovered", {
  metrics <- NULL

  expect_warning(
    metrics <- evaluate_smooth_continuous(),
    regexp = NA
  )

  expect_true(all(is.finite(unlist(metrics))))

  expect_lt(metrics$gp_normalized_rmse, 0.01)
  expect_lt(metrics$gpp_normalized_rmse, 0.01)

  expect_lt(metrics$gp_relative_rms, 0.02)
  expect_lt(metrics$gpp_relative_rms, 0.02)

  expect_lte(metrics$density_relative_l1, 0.10)

  expect_lte(
    metrics$total_modulus_relative_error,
    0.02
  )

  expect_lte(metrics$center_log_error, 0.05)

  expect_lte(
    abs(
      metrics$estimated_log_width /
        metrics$true_log_width - 1
    ),
    0.10
  )

  expect_lte(
    abs(metrics$truth_mass_on_fit_grid - 1),
    1e-3
  )
})
