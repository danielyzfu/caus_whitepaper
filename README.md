# CUPED White Paper (ECMA 31370)

## Simulation DGP

The proof-of-concept simulation models a Hulu A/B test. Full code is in `cuped_simulation.Rmd`.

### Setup

- 50,000 users, randomly assigned to treatment/control
- 20% are new users (no pre-experiment data)
- 80% are existing users with 4 weeks of pre-experiment streaming data

### How weekly streaming hours are generated

Each user has a **fixed baseline** `mu_i ~ LogNormal(log(7), 0.7)`. This gives a right-skewed distribution with median ~7 hrs/week, heavy users up to 30+.

Weekly hours for user i in week t:

```
hours_it = max(0, mu_i + tau * W_i + noise_it)
```

- `mu_i` is the same across all weeks (this is what creates pre-post correlation)
- `noise_it ~ N(0, 5^2)` is independent week-to-week noise
- `tau = 0.10` hrs/week (true treatment effect, ~1% relative lift)
- `max(0, ...)` floors at zero (no negative streaming hours)

### Why this DGP is interesting

1. **Pre-post correlation emerges naturally.** Because pre and post weeks share the same `mu_i`, they're correlated without hardcoding rho. The resulting rho is ~0.70 for a 1-week window, ~0.80 for 4-week.

2. **Zero-floor creates non-linearity.** Light users with small `mu_i` get pushed to zero by noise. These correlated zeros mean the relationship between pre and post isn't perfectly linear, so the optimal theta is not 1. It ranges from ~0.69 (1-week window) to ~0.90 (4-week window). This makes First Differences (theta=1) and CUPED (optimal theta) genuinely different without needing artificial setup.

3. **New users with missing data.** 20% of users have no pre-experiment observations. We handle them with the Booking.com COALESCE approach (leave unadjusted).

4. **Multiple covariate windows.** We construct 1-week, 2-week, and 4-week pre-experiment averages to show diminishing returns on window length (per Faire's findings).

### Methods compared

- **Naive DM**: difference-in-means, ignores pre-data
- **First Differences**: Y_post - X_pre (theta forced to 1)
- **CUPED**: optimal theta from data
- **OLS Regression Adjustment**: lm(Y ~ W + X), identical to CUPED with a single covariate

### Key results (existing users, 4-week window)

- Naive DM: cannot detect the 0.10 effect (CI includes zero)
- First Diff: significant, ~62% variance reduction
- CUPED: significant, ~63% variance reduction
- CUPED = OLS when using a single covariate

### For Monte Carlo simulations

If you want to run repeated sampling with this DGP, the `generate_week()` function in `cuped_simulation.Rmd` can be called in a loop. The DGP naturally produces method differences (Naive vs First Diff vs CUPED) without needing to add irrelevant covariates.

## File structure

- `White_Paper_Template.Rmd` - main paper (Sections 1, 2, 5, 7 drafted)
- `cuped_simulation.Rmd` - standalone simulation code (runnable)
- `whitepaper_sim.Rmd` - Monte Carlo comparison simulation
- `comparisons.png` - density plot from Monte Carlo
