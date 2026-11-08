for (seed_value in c(123L, 456L, 789L)) {
  local({
    seed <- seed_value

    test_that(
      paste("three Maxwell modes are recovered with 1% noise, seed", seed),
      {
        metrics <- NULL

        expect_warning(
          metrics <- evaluate_noisy_maxwell(
            seed = seed,
            noise_relative = 0.01
          ),
          regexp = NA
        )

        expect_equal(metrics$distinct_matches, 3L)
        expect_lt(metrics$gp_error_to_truth, 0.01)
        expect_lt(metrics$gpp_error_to_truth, 0.01)

        expect_lte(
          metrics$total_modulus_relative_error,
          0.02
        )

        expect_lte(metrics$max_log_tau_error, 0.10)

        expect_lte(
          metrics$max_mode_relative_error,
          0.10
        )
      }
    )
  })
}
