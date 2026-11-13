# ==============================================================================
# 05  COLONIZER MODEL
#     A beta binomial model of the colonizer share, fitted by maximum likelihood
# ==============================================================================
#
# THIS IS THE CUSTOM MODEL OF THE PROJECT
#   The course asks for a model-based analysis using maximum likelihood or Bayesian
#   inference. This script writes the likelihood by hand (no model-fitting package),
#   fits it with optim(), and gets confidence intervals from the profile likelihood,
#   the same approach as the Lab 6 LD50 exercise.
#
# QUESTION THIS ANSWERS
#   Is the difference between heavily impacted and unimpacted communities carried by
#   colonizer nematodes (c-p class 1), and did the colonizer share change differently
#   at heavily impacted stations between 2010 and 2011?
#   This question arose from exploring the data (script 04), so treat the results
#   as exploratory and not as a confirmatory test.
#
# THE MODEL
#   y_st ~ BetaBinomial(n_st, p_st, kappa)
#       y = number of colonizer (c-p 1) animals at station s in year t
#       n = total animals counted there (both layers added together)
#       p = the expected colonizer proportion
#       kappa = how tightly stations cluster around that proportion. Small kappa =
#               stations of the same impact level differ a lot from each other.
#               Very large kappa makes the model an ordinary binomial.
#
#   logit(p_st) = a + bHI*HI_s + bMI*MI_s + (c + dHI*HI_s + dMI*MI_s) * (year_t - 2010)
#       a    log odds of colonizers at NI (unimpacted) stations in 2010
#       bHI  how much higher (or lower) the HI log odds are than NI in 2010
#       bMI  same for MI
#       c    yearly change in log odds at NI stations
#       dHI  how much HI's yearly change differs from NI's  <-- the key parameter.
#            If HI were converging on the controls, their share would fall relative
#            to NI and dHI would be negative.
#       dMI  same for MI
#   The logit link keeps the proportion between 0 and 1 and lets us read effects as
#   odds ratios (exp of the log odds).
#
# WHY BETA BINOMIAL INSTEAD OF BINOMIAL
#   A plain binomial treats every animal as an independent draw with the same chance
#   of being a colonizer. Stations of the same impact level clearly differ by more
#   than that allows, which would make the model overconfident. The beta binomial
#   lets each station have its own true proportion. A plain binomial is also fitted,
#   and the two are compared with AIC.
#
# A DATA CAVEAT
#   2011 gives real counts. 2010 gives percentages only, so the number of animals in
#   each 2010 sample is reconstructed (see recover_n). Because that is approximate,
#   the model is refitted with all reconstructed 2010 sizes halved and doubled to see
#   whether the answer changes (it barely does).
#
# RUN
#   From the project root, after 01_data_prep.R. The profile likelihood step refits
#   the model many times, so this script takes a minute or two.
#
# INPUT
#   data/processed/community_genus_long.csv
#
# OUTPUTS
#   outputs/tables/colonizer_model_data.csv         the counts that go into the model
#   outputs/tables/colonizer_model_estimates.csv    estimates, intervals, odds ratios (Table 3)
#   outputs/tables/colonizer_model_fit.csv          beta binomial vs binomial (AIC)
#   outputs/tables/colonizer_model_sensitivity.csv  halved / doubled 2010 sample sizes
#   outputs/tables/colonizer_model_fitted.csv       fitted proportions by group and year
#   outputs/figures/colonizer_profiles.png          profile likelihoods (Figure 5)
#   outputs/figures/colonizer_fit.png               data and fitted model (Figure 4)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

genus <- read_csv("data/processed/community_genus_long.csv", show_col_types = FALSE)

# ---- 1. Build counts -------------------------------------------------------
# Reconstruct how many animals were in a 2010 sample. Percentages are counts / n x 100,
# so the right n is the one that turns every percentage back into (nearly) a whole
# number of animals. We try n = 20, 21, ... 300 and take the first that works to
# within `tol` of a whole number. If none works, fall back to assuming the rarest
# taxon was a single animal (n = 100 / smallest percentage).
recover_n <- function(p, lo = 20, hi = 300, tol = 0.1) {
  p <- p[p > 0]
  for (n in lo:hi) {
    if (max(abs(p * n / 100 - round(p * n / 100))) < tol) return(n)
  }
  round(100 / min(p))
}

# Build the 18 rows of model data (9 stations x 2 years): colonizers y out of
# total n. `scale2010` multiplies the reconstructed 2010 sample sizes (1 = as
# reconstructed; 0.5 and 2 are used in the sensitivity check).
build_data <- function(scale2010 = 1) {
  # Reconstructed size of each 2010 sample.
  n10 <- genus |>
    filter(year == 2010) |>
    summarise(n_sample = recover_n(raw_value) * scale2010, .by = sample)

  # Convert every record to a count: 2010 percentages x sample size / 100;
  # 2011 values are already counts.
  counts <- genus |>
    left_join(n10, by = "sample") |>
    mutate(count = if_else(year == 2010, raw_value * n_sample / 100, raw_value))

  # Add the two layers: y = colonizers (c-p 1), n = all animals, per station-year.
  # Indicator columns HI and MI (0/1) and year0 (0 = 2010, 1 = 2011) feed the model;
  # NI is the baseline when both indicators are 0.
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
# `par` is the parameter vector, in this order:
#   par[1] a, par[2] bHI, par[3] bMI, par[4] c, par[5] dHI, par[6] dMI, par[7] log(kappa)
# kappa is stored on the log scale so the optimiser can never make it negative.

# Linear predictor: the logit of the colonizer proportion for every row of `d`.
lin_pred <- function(par, d) {
  par[1] + par[2] * d$HI + par[3] * d$MI +
    (par[4] + par[5] * d$HI + par[6] * d$MI) * d$year0
}

# Negative log likelihood of the beta binomial model (what optim() minimises).
# Each station's true proportion is drawn from Beta(p*kappa, (1-p)*kappa), then y
# is binomial given that proportion. Integrating that out gives, per observation,
#   log choose(n, y) + log B(y + a, n - y + b) - log B(a, b)
# with a = p*kappa and b = (1-p)*kappa, where B is the beta function. We use
# lchoose() and lbeta() (the log versions) for numerical stability.
nll_bb <- function(par, d) {
  p <- plogis(lin_pred(par, d))     # inverse logit: back to a proportion
  kappa <- exp(par[7])
  a <- p * kappa
  b <- (1 - p) * kappa
  -sum(lchoose(d$n, d$y) + lbeta(d$y + a, d$n - d$y + b) - lbeta(a, b))
}

# Negative log likelihood of an ordinary binomial model (no extra variation),
# fitted only so we can compare it with the beta binomial. Uses 6 parameters.
nll_binom <- function(par, d) {
  -sum(dbinom(d$y, d$n, plogis(lin_pred(par, d)), log = TRUE))
}

# Starting values for the optimiser: `a` starts at the observed NI colonizer
# proportion in 2010 (the +0.5 and +1 avoid taking the log of zero); every effect
# starts at 0 (no difference); kappa starts at 20 (moderate extra variation).
start_for <- function(d) {
  ni <- d$impact == "NI" & d$year0 == 0
  c(a = qlogis((sum(d$y[ni]) + 0.5) / (sum(d$n[ni]) + 1)), bHI = 0, bMI = 0, c = 0, dHI = 0, dMI = 0,
    logkappa = log(20))
}

# Minimise the negative log likelihood with BFGS (a quasi-Newton method).
# hessian = TRUE also returns the curvature at the minimum, used below for standard
# errors and for the width of the profile grid.
fit_bb <- function(d) {
  optim(start_for(d), nll_bb, d = d, method = "BFGS", hessian = TRUE, control = list(maxit = 2000))
}

fit <- fit_bb(dat)
# convergence == 0 means optim() reports success; stop if it did not.
stopifnot(fit$convergence == 0)

# The binomial comparison model (6 parameters: no kappa).
fit_binom <- optim(start_for(dat)[1:6], nll_binom, d = dat, method = "BFGS", control = list(maxit = 2000))

# AIC = 2 x (negative log likelihood) + 2 x (number of parameters). It rewards
# fit but penalises extra parameters; lower is better.
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

# Wald standard errors: square roots of the diagonal of the inverse Hessian. These
# are quick but assume the likelihood is symmetric near its peak, so we use them
# only to size the profile grid and as a comparison; the reported intervals are
# profile likelihood intervals.
se <- sqrt(diag(solve(fit$hessian)))

# The 95% profile cutoff. A parameter value is inside the 95% interval if the
# negative log likelihood, re-optimised over the other parameters, is no more than
# qchisq(0.95, 1) / 2 = 1.92 above the minimum.
cutoff_drop <- qchisq(0.95, df = 1) / 2

# PROFILE LIKELIHOOD for one parameter (`idx`): fix it at each of 61 values spread
# over +/- 5 standard errors of the estimate, and at each value re-fit all the OTHER
# parameters. The resulting curve shows how fast the fit worsens as that parameter is
# pushed away from its best value; it need not be symmetric.
profile_par <- function(idx, d = dat, mle = fit$par, grid_se = fit$se_use) {
  grid <- seq(mle[idx] - 5 * grid_se[idx], mle[idx] + 5 * grid_se[idx], length.out = 61)
  nll <- vapply(grid, function(v) {
    # Objective over the remaining 6 parameters, with parameter `idx` held at v.
    f <- function(rest) {
      par <- numeric(7)
      par[idx] <- v
      par[-idx] <- rest
      nll_bb(par, d)
    }
    # Start the search from the overall best values (minus the fixed parameter).
    optim(mle[-idx], f, method = "BFGS", control = list(maxit = 1000))$value
  }, numeric(1))
  tibble(idx = idx, parameter = par_names[idx], value = grid, nll = nll)
}

# (profile_par reads the standard errors from fit$se_use, so store them there.)
fit$se_use <- se
# Profile the six regression parameters (log kappa is not of interest).
profiles <- bind_rows(lapply(1:6, profile_par))

# Read the 95% interval off a profile curve: find where the curve crosses the cutoff
# on each side of the best estimate, by straight line interpolation between grid
# points. NA means the curve never crossed within the grid (interval would extend
# beyond the grid).
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

# One row per parameter: estimate, Wald standard error and interval, and the
# profile interval; then the odds ratio version (exp of the log odds).
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
# Refit the whole model with every reconstructed 2010 sample size halved, as
# reconstructed, and doubled. If bHI and dHI hardly move, the conclusions do not
# depend on the exact reconstruction.
sens <- lapply(c(0.5, 1, 2), function(s) {
  f <- fit_bb(build_data(s))
  tibble(scale_2010_n = s, converged = f$convergence == 0,
         bHI = f$par[2], dHI = f$par[5], kappa = exp(f$par[7]))
}) |> bind_rows()
write_csv(sens, "outputs/tables/colonizer_model_sensitivity.csv")

# ---- 5. Figures ------------------------------------------------------------
# Figure 5: profile curves for the two parameters the paper focuses on (bHI and
# dHI). The 95% interval is where the curve is below the dashed cutoff line; the
# grey vertical line marks the best estimate. The y axis is capped at 6 so the
# region around the cutoff is easy to see.
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

# Figure 4: fitted colonizer proportion for each impact group in each year (the
# model's lines), with the observed station proportions as points.
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
