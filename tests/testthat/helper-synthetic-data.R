make_maxwell_data <- function(
	omega,
	G = c(1500, 400, 80),
	tau = c(0.005, 0.1, 2),
	noise_relative = 0,
	seed = 123
) {
	stopifnot(
		is.numeric(omega),
		all(is.finite(omega)),
		all(omega > 0),
		is.numeric(G),
		is.numeric(tau),
		length(G) == length(tau),
		all(is.finite(G)),
		all(is.finite(tau)),
		all(G > 0),
		all(tau > 0),
		length(noise_relative) == 1L,
		is.finite(noise_relative),
		noise_relative >= 0
	)

	wt <- outer(omega, tau, "*")

	Gp <- rowSums(
		sweep(
			wt^2 / (1 + wt^2),
			2,
			G,
			"*"
		)
	)

	Gpp <- rowSums(
		sweep(
			wt / (1 + wt^2),
			2,
			G,
			"*"
		)
	)

	if (noise_relative > 0) {
		set.seed(seed)

		Gp <- Gp * (
			1 + rnorm(
				length(Gp),
				mean = 0,
				sd = noise_relative
			)
		)

		Gpp <- Gpp * (
			1 + rnorm(
				length(Gpp),
				mean = 0,
				sd = noise_relative
			)
		)
	}

	data.frame(
		omega = omega,
		Gp = Gp,
		Gpp = Gpp
	)
}
