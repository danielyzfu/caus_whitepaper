# =============================================================================
# CUPED Simulation: Streaming Platform A/B Test
# =============================================================================
# Simulates a mid-tier streaming platform running an A/B test on a new
# recommendation algorithm.
#
# Data generating process:
#   1. Each user has a stable baseline streaming rate (lognormal, right-skewed)
#   2. Weekly hours = max(0, baseline + noise)
#      Light users naturally produce zero-hours weeks when noise pushes below 0.
#      These zeros are correlated across weeks (same users tend to be inactive),
#      preserving realistic pre-post correlation.
#   3. Treatment group gets a small boost in streaming hours
#   4. New users (20%) have no pre-experiment data
#
# Target properties (aligned with Netflix findings):
#   - Pre-post rho ≈ 0.70 for 1-week window (existing users)
#   - Variance reduction ≈ 35-45% for existing users (matches Netflix ~40%)
#   - Naive method: borderline / not significant
#   - CUPED: significant
#   - Longer pre-experiment window → higher rho → more reduction
# =============================================================================

library(tidyverse)

set.seed(42)

# =============================================================================
# 1. Parameters
# =============================================================================

N <- 50000
new_user_frac <- 0.20

# baseline streaming hours (lognormal)
log_mean <- log(7)       # median ~7 hrs/week
log_sd <- 0.7            # right-skewed range: ~1 to 30+ hrs/week

# week-to-week noise (larger → lower rho, more zeros for light users)
noise_sd <- 5.0

# treatment effect: ~1.5% relative lift on mean of ~9 hrs/week
# small enough that naive method struggles, large enough that CUPED catches it
tau <- 0.10

n_pre_weeks <- 4

# =============================================================================
# 2. Generate Users
# =============================================================================

n_new <- round(N * new_user_frac)
n_existing <- N - n_new

users <- tibble(
  user_id = 1:N,
  is_new = c(rep(TRUE, n_new), rep(FALSE, n_existing)),
  baseline_mu = rlnorm(N, meanlog = log_mean, sdlog = log_sd),
  treatment = sample(c(0L, 1L), N, replace = TRUE)
)

# =============================================================================
# 3. Generate Weekly Streaming Hours
# =============================================================================

generate_week <- function(baseline_mu, treatment_indicator = NULL, tau_val = 0) {
  n <- length(baseline_mu)
  noise <- rnorm(n, mean = 0, sd = noise_sd)
  hours <- baseline_mu + noise
  if (!is.null(treatment_indicator)) {
    hours <- hours + tau_val * treatment_indicator
  }
  return(pmax(hours, 0))
}

pre_weeks <- map(1:n_pre_weeks, ~ generate_week(users$baseline_mu))
names(pre_weeks) <- paste0("pre_week_", 1:n_pre_weeks)

post_week <- generate_week(users$baseline_mu, users$treatment, tau)

# =============================================================================
# 4. Assemble Dataset
# =============================================================================

sim_data <- users %>%
  bind_cols(as_tibble(pre_weeks)) %>%
  mutate(
    post_week = post_week,
    pre_1wk = pre_week_4,
    pre_2wk = (pre_week_3 + pre_week_4) / 2,
    pre_4wk = (pre_week_1 + pre_week_2 + pre_week_3 + pre_week_4) / 4
  ) %>%
  mutate(
    pre_1wk = if_else(is_new, NA_real_, pre_1wk),
    pre_2wk = if_else(is_new, NA_real_, pre_2wk),
    pre_4wk = if_else(is_new, NA_real_, pre_4wk)
  )

# =============================================================================
# 5. Data Inspection
# =============================================================================

existing <- sim_data %>% filter(!is_new)

cat("===========================================================\n")
cat("  DATA GENERATION SUMMARY\n")
cat("===========================================================\n\n")

cat(sprintf("Total users: %d\n", nrow(sim_data)))
cat(sprintf("  Treatment: %d | Control: %d\n",
            sum(sim_data$treatment == 1), sum(sim_data$treatment == 0)))
cat(sprintf("  New: %d (%.0f%%) | Existing: %d (%.0f%%)\n",
            n_new, new_user_frac * 100, n_existing, (1 - new_user_frac) * 100))
cat(sprintf("  True tau = %.2f hrs/week\n\n", tau))

cat("--- Streaming Hours (post-experiment) ---\n")
cat(sprintf("  Mean:   %.2f hrs/week\n", mean(sim_data$post_week)))
cat(sprintf("  Median: %.2f hrs/week\n", median(sim_data$post_week)))
cat(sprintf("  SD:     %.2f\n", sd(sim_data$post_week)))
pct_zero <- mean(sim_data$post_week == 0) * 100
cat(sprintf("  Zero-hours: %.1f%% of users\n\n", pct_zero))

rho_1 <- cor(existing$pre_1wk, existing$post_week)
rho_2 <- cor(existing$pre_2wk, existing$post_week)
rho_4 <- cor(existing$pre_4wk, existing$post_week)

cat("--- Pre-Post Correlation (existing users) ---\n")
cat(sprintf("  1-week:  rho = %.3f  (var reduction ≈ %.1f%%)\n", rho_1, rho_1^2 * 100))
cat(sprintf("  2-week:  rho = %.3f  (var reduction ≈ %.1f%%)\n", rho_2, rho_2^2 * 100))
cat(sprintf("  4-week:  rho = %.3f  (var reduction ≈ %.1f%%)\n", rho_4, rho_4^2 * 100))
cat("\n")

# =============================================================================
# 6. CUPED (Booking.com COALESCE approach for missing data)
# =============================================================================

cuped_adjust <- function(y, x) {
  has_x <- !is.na(x)
  theta <- cov(y[has_x], x[has_x]) / var(x[has_x])
  x_mean <- mean(x[has_x])
  y_adj <- y
  y_adj[has_x] <- y[has_x] - theta * (x[has_x] - x_mean)
  return(list(y_adj = y_adj, theta = theta))
}

# =============================================================================
# 7. Naive vs CUPED
# =============================================================================

run_test <- function(data, pre_col, label) {
  y <- data$post_week
  x <- data[[pre_col]]
  w <- data$treatment

  # naive
  y_t <- y[w == 1]; y_c <- y[w == 0]
  n_t <- length(y_t); n_c <- length(y_c)
  naive_est <- mean(y_t) - mean(y_c)
  naive_se <- sqrt(var(y_t) / n_t + var(y_c) / n_c)
  naive_ci <- naive_est + c(-1, 1) * 1.96 * naive_se
  naive_sig <- (naive_ci[1] > 0) | (naive_ci[2] < 0)

  # cuped
  res <- cuped_adjust(y, x)
  ya <- res$y_adj
  ya_t <- ya[w == 1]; ya_c <- ya[w == 0]
  cuped_est <- mean(ya_t) - mean(ya_c)
  cuped_se <- sqrt(var(ya_t) / n_t + var(ya_c) / n_c)
  cuped_ci <- cuped_est + c(-1, 1) * 1.96 * cuped_se
  cuped_sig <- (cuped_ci[1] > 0) | (cuped_ci[2] < 0)

  vr <- 1 - cuped_se^2 / naive_se^2

  cat(sprintf("  [%s]\n", label))
  cat(sprintf("    Naive: est=%.4f  SE=%.4f  CI=[%.4f, %.4f]  %s\n",
              naive_est, naive_se, naive_ci[1], naive_ci[2],
              ifelse(naive_sig, "SIG", "n.s.")))
  cat(sprintf("    CUPED: est=%.4f  SE=%.4f  CI=[%.4f, %.4f]  %s\n",
              cuped_est, cuped_se, cuped_ci[1], cuped_ci[2],
              ifelse(cuped_sig, "SIG", "n.s.")))
  cat(sprintf("    Var reduction: %.1f%%\n\n", vr * 100))

  return(tibble(label, naive_est, naive_se, naive_sig, cuped_est, cuped_se, cuped_sig, vr))
}

cat("===========================================================\n")
cat("  RESULTS: ALL USERS\n")
cat("===========================================================\n\n")

r_all <- bind_rows(
  run_test(sim_data, "pre_1wk", "All users, 1-wk window"),
  run_test(sim_data, "pre_2wk", "All users, 2-wk window"),
  run_test(sim_data, "pre_4wk", "All users, 4-wk window")
)

cat("===========================================================\n")
cat("  RESULTS: EXISTING USERS ONLY\n")
cat("===========================================================\n\n")

r_ex <- bind_rows(
  run_test(existing, "pre_1wk", "Existing, 1-wk window"),
  run_test(existing, "pre_2wk", "Existing, 2-wk window"),
  run_test(existing, "pre_4wk", "Existing, 4-wk window")
)

cat("===========================================================\n")
cat("  SUBGROUP: NEW USERS\n")
cat("===========================================================\n\n")

new_u <- sim_data %>% filter(is_new)
yt <- new_u$post_week[new_u$treatment == 1]
yc <- new_u$post_week[new_u$treatment == 0]
est <- mean(yt) - mean(yc)
se <- sqrt(var(yt) / length(yt) + var(yc) / length(yc))
ci <- est + c(-1, 1) * 1.96 * se
sig <- (ci[1] > 0) | (ci[2] < 0)
cat(sprintf("  Naive: est=%.4f  SE=%.4f  CI=[%.4f, %.4f]  %s\n",
            est, se, ci[1], ci[2], ifelse(sig, "SIG", "n.s.")))
cat("  (No CUPED possible for new users)\n\n")

# =============================================================================
# 8. Save
# =============================================================================

output <- sim_data %>%
  select(user_id, is_new, treatment, pre_1wk, pre_2wk, pre_4wk, post_week)

write_csv(output, "/tmp/caus_whitepaper/code/simulated_streaming_data.csv")
cat(sprintf("Saved: %d rows x %d cols → simulated_streaming_data.csv\n",
            nrow(output), ncol(output)))
