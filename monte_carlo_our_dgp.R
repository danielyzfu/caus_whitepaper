library(tidyverse)
set.seed(54321)

# Our Hulu DGP parameters
log_mean <- log(7)
log_sd <- 0.7
noise_sd <- 5.0
tau <- 0.10
n <- 500  # moderate sample size per run

r <- 1000  # number of Monte Carlo repetitions

estimates <- map_dfr(1:r, function(i) {
  # Generate users
  baseline_mu <- rlnorm(n, meanlog = log_mean, sdlog = log_sd)
  w <- rbinom(n, 1, 0.5)

  # Pre-experiment week (no treatment)
  pre <- pmax(0, baseline_mu + rnorm(n, 0, noise_sd))

  # Post-experiment week (with treatment)
  post <- pmax(0, baseline_mu + tau * w + rnorm(n, 0, noise_sd))

  # Naive DM
  naive_est <- mean(post[w == 1]) - mean(post[w == 0])

  # First Differences (theta = 1)
  fd <- post - pre
  fd_est <- mean(fd[w == 1]) - mean(fd[w == 0])

  # CUPED (optimal theta)
  theta <- cov(post, pre) / var(pre)
  y_adj <- post - theta * (pre - mean(pre))
  cuped_est <- mean(y_adj[w == 1]) - mean(y_adj[w == 0])

  data.frame(
    Naive = naive_est,
    FirstDiff = fd_est,
    CUPED = cuped_est
  )
})

# Plot
estimates_long <- estimates %>%
  pivot_longer(everything(), names_to = "Estimator", values_to = "Estimate")

estimates_long$Estimator <- factor(estimates_long$Estimator,
                                    levels = c("Naive", "FirstDiff", "CUPED"))

p <- ggplot(estimates_long, aes(x = Estimate, color = Estimator)) +
  geom_density(linewidth = 0.8) +
  geom_vline(xintercept = tau, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Naive" = "gray50", "FirstDiff" = "orange3", "CUPED" = "steelblue")) +
  labs(title = "Monte Carlo Sampling Distributions (Hulu DGP, n=500, 1000 reps)",
       x = "Estimated Treatment Effect",
       y = "Density") +
  theme_minimal()

ggsave("/tmp/caus_whitepaper/code/mc_hulu_dgp.png", plot = p, width = 8, height = 5, dpi = 300)
cat("Done. Saved to mc_hulu_dgp.png\n")

# Print SDs
cat(sprintf("\nSD of estimates:\n"))
cat(sprintf("  Naive:     %.4f\n", sd(estimates$Naive)))
cat(sprintf("  FirstDiff: %.4f\n", sd(estimates$FirstDiff)))
cat(sprintf("  CUPED:     %.4f\n", sd(estimates$CUPED)))
