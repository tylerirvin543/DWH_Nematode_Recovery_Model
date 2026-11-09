# 06_colonizer_model.R
# A hand-written beta-binomial model for the proportion of colonizer (c-p 1) nematodes,
# fitted by maximum likelihood with optim, with profile likelihood intervals.
#
# Run from the project root, after scripts/01_data_prep.R.
#
# Inputs : data/processed/community_genus_long.csv
# Outputs: outputs/tables/colonizer_model_data.csv, colonizer_model_estimates.csv,
#          colonizer_model_fit.csv, colonizer_model_sensitivity.csv, colonizer_model_fitted.csv
#          outputs/figures/colonizer_profiles.png, colonizer_fit.png
#
# Model
#   y_st ~ BetaBinomial(n_st, p_st, kappa)       y = c-p 1 animals, n = all animals (station-year)
#   logit(p_st) = a + bHI*HI_s + bMI*MI_s + (c + dHI*HI_s + dMI*MI_s) * (year_t - 2010)
#   a    : log odds of colonizers at unimpacted (NI) stations in 2010
#   bHI  : HI minus NI difference in log odds in 2010 (bMI likewise)
#   c    : yearly change in log odds at NI stations
#   dHI  : how much more the HI yearly change differs from NI's (dMI likewise)
#   kappa: concentration; small kappa means extra variation between stations
#
# 2011 gives counts. 2010 gives percentages only, so 2010 sample sizes are reconstructed
# (see recover_n) and the model is refitted with those sizes halved and doubled as a check.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

genus <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE)

# ---- 1. Build counts -------------------------------------------------------
# Find the sample size that turns a sample's percentages into whole numbers.
# Falls back to 100 / (smallest percentage), i.e. assuming the rarest taxon is a singleton.
recover_n <- function(p, lo = 20, hi = 300, tol = 0.1) {
  p <- p[p > 0]
  for (n in lo:hi) {
    if (max(abs(p * n / 100 - round(p * n / 100))) < tol) return(n)
  }
  round(100 / min(p))
}

build_data <- function(scale2010 = 1) {
  n10 <- genus |>
    filter(year == 2010) |>
    summarise(n_sample = recover_n(raw_value) * scale2010, .by = sample)

  counts <- genus |>
    left_join(n10, by = "sample") |>
    mutate(count = if_else(year == 2010, raw_value * n_sample / 100, raw_value))

  counts |>
    summarise(
      y = round(sum(count[cp == 1])),
      n = round(sum(count)),
      .by = c(year, station, impact)
    ) |>
    mutate(
      HI = as.numeric(impact == "HI"),
      MI = as.numeric(impact == "MI"),
      year0 = year - 2010,
      impact = factor(impact, levels = c("HI", "MI", "NI"))
    ) |>
    arrange(year, impact, station)
}

dat <- build_data()
stopifnot(nrow(dat) == 18)
write_csv(dat |> select(year, station, impact, y, n), "outputs/tables/colonizer_model_data.csv")

# ---- 2. Likelihoods --------------------------------------------------------
lin_pred <- function(par, d) {
  par[1] + par[2] * d$HI + par[3] * d$MI +
    (par[4] + par[5] * d$HI + par[6] * d$MI) * d$year0
}

# Beta-binomial: p ~ Beta(p*kappa, (1-p)*kappa), then y ~ Binomial(n, p)
nll_bb <- function(par, d) {
  p <- plogis(lin_pred(par, d))
  kappa <- exp(par[7])
  a <- p * kappa
  b <- (1 - p) * kappa
  -sum(lchoose(d$n, d$y) + lbeta(d$y + a, d$n - d$y + b) - lbeta(a, b))
}

# Plain binomial (no extra variation), for comparison
nll_binom <- function(par, d) {
  -sum(dbinom(d$y, d$n, plogis(lin_pred(par, d)), log = TRUE))
}

start_for <- function(d) {
  ni <- d$impact == "NI" & d$year0 == 0
  c(a = qlogis((sum(d$y[ni]) + 0.5) / (sum(d$n[ni]) + 1)), bHI = 0, bMI = 0, c = 0, dHI = 0, dMI = 0,
    logkappa = log(20))
}

fit_bb <- function(d) {
  optim(start_for(d), nll_bb, d = d, method = "BFGS", hessian = TRUE, control = list(maxit = 2000))
}

fit <- fit_bb(dat)
stopifnot(fit$convergence == 0)

fit_binom <- optim(start_for(dat)[1:6], nll_binom, d = dat, method = "BFGS", control = list(maxit = 2000))

fit_table <- tibble(
  model = c("Beta-binomial", "Binomial"),
  n_parameters = c(7, 6),
  neg_log_lik = c(fit$value, fit_binom$value),
  AIC = 2 * c(fit$value, fit_binom$value) + 2 * c(7, 6),
  converged = c(fit$convergence == 0, fit_binom$convergence == 0)
)
write_csv(fit_table, "outputs/tables/colonizer_model_fit.csv")

# ---- 3. Intervals ----------------------------------------------------------
par_names <- c("a (NI log odds, 2010)", "bHI (HI minus NI, 2010)", "bMI (MI minus NI, 2010)",
               "c (NI change per year)", "dHI (HI minus NI, change per year)",
               "dMI (MI minus NI, change per year)", "log kappa")

se <- sqrt(diag(solve(fit$hessian)))
cutoff_drop <- qchisq(0.95, df = 1) / 2   # 1.92 log-likelihood units

# Profile likelihood: fix one parameter on a grid, re-optimise the others
profile_par <- function(idx, d = dat, mle = fit$par, grid_se = fit$se_use) {
  grid <- seq(mle[idx] - 5 * grid_se[idx], mle[idx] + 5 * grid_se[idx], length.out = 61)
  nll <- vapply(grid, function(v) {
    f <- function(rest) {
      par <- numeric(7)
      par[idx] <- v
      par[-idx] <- rest
      nll_bb(par, d)
    }
    optim(mle[-idx], f, method = "BFGS", control = list(maxit = 1000))$value
  }, numeric(1))
  tibble(idx = idx, parameter = par_names[idx], value = grid, nll = nll)
}

fit$se_use <- se
profiles <- bind_rows(lapply(1:6, profile_par))

profile_ci <- function(pr, mle_val, min_nll) {
  lim <- min_nll + cutoff_drop
  edge <- function(side) {
    s <- if (side == "lo") pr$value <= mle_val else pr$value >= mle_val
    x <- pr$nll[s]; y <- pr$value[s]
    if (!any(x > lim)) return(NA_real_)
    approx(x, y, xout = lim, ties = mean)$y
  }
  c(lo = edge("lo"), hi = edge("hi"))
}

estimates <- lapply(1:6, function(i) {
  pr <- profiles |> filter(idx == i)
  ci <- profile_ci(pr, fit$par[i], fit$value)
  est <- unname(fit$par[i])
  s <- unname(se[i])
  tibble(
    parameter = par_names[i],
    estimate = est, std_error = s,
    wald_lower = est - 1.96 * s, wald_upper = est + 1.96 * s,
    profile_lower = unname(ci["lo"]), profile_upper = unname(ci["hi"])
  )
}) |> bind_rows() |>
  mutate(odds_ratio = exp(estimate),
         or_lower = exp(profile_lower), or_upper = exp(profile_upper))

write_csv(estimates, "outputs/tables/colonizer_model_estimates.csv")

# ---- 4. Sensitivity to the reconstructed 2010 sample sizes -----------------
sens <- lapply(c(0.5, 1, 2), function(s) {
  f <- fit_bb(build_data(s))
  tibble(scale_2010_n = s, converged = f$convergence == 0,
         bHI = f$par[2], dHI = f$par[5], kappa = exp(f$par[7]))
}) |> bind_rows()
write_csv(sens, "outputs/tables/colonizer_model_sensitivity.csv")

# ---- 5. Figures ------------------------------------------------------------
prof_plot_df <- profiles |>
  filter(idx %in% c(2, 5)) |>
  mutate(parameter = factor(parameter, levels = par_names[c(2, 5)]))

mle_df <- tibble(parameter = factor(par_names[c(2, 5)], levels = par_names[c(2, 5)]),
                 est = fit$par[c(2, 5)])

prof_plot <- ggplot(prof_plot_df, aes(value, nll - fit$value)) +
  geom_line() +
  geom_hline(yintercept = cutoff_drop, linetype = "dashed") +
  geom_vline(data = mle_df, aes(xintercept = est), colour = "grey50") +
  facet_wrap(~parameter, scales = "free_x") +
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = "Parameter value (log odds)", y = "Increase in negative log likelihood",
       subtitle = "Dashed line marks the 95% profile likelihood cutoff") +
  theme_bw()

ggsave("outputs/figures/colonizer_profiles.png", prof_plot, width = 8, height = 4, dpi = 300)

fitted_df <- expand_grid(impact = factor(c("HI", "MI", "NI"), levels = c("HI", "MI", "NI")),
                         year = c(2010, 2011)) |>
  mutate(HI = as.numeric(impact == "HI"), MI = as.numeric(impact == "MI"), year0 = year - 2010)
fitted_df$p <- plogis(lin_pred(fit$par, fitted_df))
write_csv(fitted_df |> select(impact, year, p), "outputs/tables/colonizer_model_fitted.csv")

fit_plot <- ggplot(dat, aes(year, y / n)) +
  geom_line(data = fitted_df, aes(y = p, colour = impact), linewidth = 1) +
  geom_point(aes(colour = impact), size = 2.5) +
  scale_x_continuous(breaks = c(2010, 2011)) +
  labs(x = "Year", y = "Proportion of colonizer (c-p 1) animals", colour = "Impact",
       subtitle = "Points are stations; lines are fitted beta-binomial means") +
  theme_bw()

ggsave("outputs/figures/colonizer_fit.png", fit_plot, width = 6.5, height = 4.5, dpi = 300)

message("Done. Colonizer model tables and figures written to outputs/.")
