###############################################################################
#  02_ols_vs_hmc.R
#  OLS regression vs. hierarchical Bayesian regression fit with HMC
#
#  Model (station j = 1..3, day i):
#    log PM2.5_ij = alpha_j + b1*WIND_ij + b2*RAIN_ij + b3*PRES_ij + e_ij
#    e_ij    ~ N(0, sigma^2)
#    alpha_j ~ N(mu_alpha, tau^2)     <- station baselines drawn from a
#                                        population distribution
#  Priors:
#    mu_alpha ~ N(0, 10^2),  b_k ~ N(0, 1),
#    tau ~ Half-Normal(0, 1),  sigma ~ Half-Normal(0, 1)
#
#  Run 01_preprocess.R first.
###############################################################################

set.seed(197)
library(coda)                      

d <- read.csv("data/air_daily.csv")
d$station <- factor(d$station)

y   <- d$log_pm25
X   <- as.matrix(d[, c("WSPM_z", "RAIN_z", "PRES_z")])
sid <- d$station_id
N <- length(y); J <- max(sid); P <- ncol(X)
stations <- levels(d$station)
pred_names <- c("Wind", "Rain", "Pressure")

# =============================================================================
# PART 1: OLS
# =============================================================================
# (a) Pooled OLS: one intercept for everyone, ignores which station a day is from
ols_pooled <- lm(log_pm25 ~ WSPM_z + RAIN_z + PRES_z, data = d)
summary(ols_pooled)

# (b) OLS with station dummies ("fixed effects"): separate intercepts, but each
#     is estimated using only that station's data, with no sharing of information.
ols_fe <- lm(log_pm25 ~ 0 + station + WSPM_z + RAIN_z + PRES_z, data = d)
summary(ols_fe)

# =============================================================================
# PART 2: HIERARCHICAL MODEL VIA HAMILTONIAN MONTE CARLO
# =============================================================================
# HMC treats the parameters as the position of a particle sliding over the
# surface U(theta) = -log posterior. It gives the particle a random momentum,
# simulates its motion (leapfrog steps), and accepts the end point with a
# Metropolis correction. It needs: (1) log posterior, (2) its gradient.
#
# Positive parameters (tau, sigma) are sampled on the log scale so that HMC
# can move freely over (-Inf, Inf).
#
# theta = [ alpha_1..alpha_J | mu_alpha | log_tau | b_1..b_P | log_sigma ]
i_alpha <- 1:J
i_mu    <- J + 1
i_ltau  <- J + 2
i_beta  <- J + 2 + 1:P
i_lsig  <- J + P + 3
D <- i_lsig
par_names <- c(paste0("alpha_", stations), "mu_alpha", "tau",
               paste0("beta_", pred_names), "sigma")

# ---- 2.1 Log posterior ------------------------------------------------------
log_post <- function(th) {
  alpha <- th[i_alpha]; mu <- th[i_mu]; beta <- th[i_beta]
  tau <- exp(th[i_ltau]); sigma <- exp(th[i_lsig])

  r <- y - alpha[sid] - X %*% beta
  ll <- -N * log(sigma) - sum(r^2) / (2 * sigma^2)                  # likelihood
  lp_alpha <- -J * log(tau) - sum((alpha - mu)^2) / (2 * tau^2)     # alpha_j ~ N(mu, tau)
  lp_prior <- -mu^2 / (2 * 100) - sum(beta^2) / 2 -                 # mu, beta priors
              tau^2 / 2 - sigma^2 / 2 +                             # half-normal(0,1)
              th[i_ltau] + th[i_lsig]                               # log-Jacobians
  ll + lp_alpha + lp_prior
}

# ---- 2.2 Gradient (derived by hand; checked numerically below) -------------
grad_log_post <- function(th) {
  alpha <- th[i_alpha]; mu <- th[i_mu]; beta <- th[i_beta]
  tau <- exp(th[i_ltau]); sigma <- exp(th[i_lsig])

  r <- as.vector(y - alpha[sid] - X %*% beta)
  dev <- alpha - mu
  g <- numeric(D)
  g[i_alpha] <- rowsum(r, sid)[, 1] / sigma^2 - dev / tau^2
  g[i_mu]    <- sum(dev) / tau^2 - mu / 100
  g[i_ltau]  <- -J + sum(dev^2) / tau^2 - tau^2 + 1
  g[i_beta]  <- crossprod(X, r) / sigma^2 - beta
  g[i_lsig]  <- -N + sum(r^2) / sigma^2 - sigma^2 + 1
  g
}

# Sanity check: analytic gradient vs. finite differences (should match ~1e-4)
th_test <- rnorm(D, 0, 0.5)
num_grad <- sapply(1:D, function(k) {
  h <- 1e-5; e <- replace(numeric(D), k, h)
  (log_post(th_test + e) - log_post(th_test - e)) / (2 * h)
})
cat("Max gradient error:", max(abs(num_grad - grad_log_post(th_test)) / (1 + abs(num_grad))), "\n")

# ---- 2.3 One HMC iteration --------------------------------------------------
# eps = leapfrog step size, L = number of steps, m_inv = inverse mass (diag),
# which rescales each parameter so they all move at a sensible speed.
hmc_step <- function(th, eps, L, m_inv) {
  p0 <- rnorm(D, 0, sqrt(1 / m_inv))                     # 1. random momentum
  th_new <- th
  p <- p0 + eps / 2 * grad_log_post(th_new)              # 2. leapfrog
  for (l in 1:L) {
    th_new <- th_new + eps * m_inv * p
    if (l < L) p <- p + eps * grad_log_post(th_new)
  }
  p <- p + eps / 2 * grad_log_post(th_new)
  H0 <- -log_post(th)     + sum(m_inv * p0^2) / 2        # 3. energy before/after
  H1 <- -log_post(th_new) + sum(m_inv * p^2)  / 2
  if (is.finite(H1) && log(runif(1)) < H0 - H1)          # 4. Metropolis accept
    list(th = th_new, acc = TRUE) else list(th = th, acc = FALSE)
}

# ---- 2.4 Full sampler with warm-up -----------------------------------------
# Warm-up (discarded, like burn-in): adapt eps toward ~65-80% acceptance and
# estimate each parameter's scale for m_inv.
run_hmc <- function(n_warm = 1500, n_samp = 2000, L = 15) {
  th <- c(rnorm(J, mean(y), 0.5), mean(y), log(0.5), rnorm(P, 0, 0.1), log(1))
  eps <- 0.005; m_inv <- rep(1, D)
  warm <- matrix(NA, n_warm, D)
  for (t in 1:n_warm) {
    s <- hmc_step(th, eps, L, m_inv); th <- s$th; warm[t, ] <- th
    eps <- eps * if (s$acc) 1.02 else 0.95
    if (t == floor(n_warm / 2)) {                  # learn scales, restart eps
      m_inv <- apply(warm[floor(n_warm / 4):t, ], 2, var) + 1e-8
      eps <- 0.1
    }
  }
  draws <- matrix(NA, n_samp, D); acc <- 0
  for (t in 1:n_samp) {
    s <- hmc_step(th, eps * runif(1, 0.8, 1.2), L, m_inv)   # jittered eps
    th <- s$th; acc <- acc + s$acc; draws[t, ] <- th
  }
  draws[, c(i_ltau, i_lsig)] <- exp(draws[, c(i_ltau, i_lsig)])  # back to tau, sigma
  colnames(draws) <- par_names
  attr(draws, "accept") <- acc / n_samp
  draws
}

chains <- lapply(1:4, function(k) run_hmc())        # 4 independent chains
cat("Acceptance rates:", round(sapply(chains, attr, "accept"), 2), "\n")
post <- mcmc.list(lapply(chains, mcmc))
draws <- as.matrix(post)

# =============================================================================
# PART 3: DIAGNOSTICS
# =============================================================================
diag_tab <- data.frame(
  Rhat = gelman.diag(post, multivariate = FALSE)$psrf[, 1],   # want < 1.01
  ESS  = effectiveSize(post))                                 # want > ~400
round(diag_tab, 3)

dir.create("figures", showWarnings = FALSE)
png("figures/hmc_traceplots.png", 1400, 900, res = 120)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
cols4 <- c("#2a78b5", "#d1603d", "#3a9d5d", "#8a5fbf")
for (v in par_names) {
  plot(NULL, xlim = c(1, nrow(chains[[1]])), ylim = range(sapply(chains, function(c) range(c[, v]))),
       main = v, xlab = "iteration", ylab = "")
  for (k in 1:4) lines(chains[[k]][, v], col = adjustcolor(cols4[k], .6))
}
dev.off()

# =============================================================================
# PART 4: RESULTS AND COMPARISON WITH OLS
# =============================================================================
summ <- t(apply(draws, 2, function(v)
  c(mean = mean(v), sd = sd(v), lo95 = quantile(v, .025, names = FALSE),
    hi95 = quantile(v, .975, names = FALSE))))
round(summ, 3)

ols_ci <- function(fit) cbind(coef(fit), confint(fit))
cmp <- rbind(
  data.frame(param = pred_names, model = "OLS pooled",
             setNames(as.data.frame(ols_ci(ols_pooled)[-1, ]), c("est", "lo", "hi"))),
  data.frame(param = pred_names, model = "OLS + station dummies",
             setNames(as.data.frame(ols_ci(ols_fe)[-(1:J), ]), c("est", "lo", "hi"))),
  data.frame(param = pred_names, model = "Hierarchical (HMC)",
             est = summ[paste0("beta_", pred_names), "mean"],
             lo  = summ[paste0("beta_", pred_names), "lo95"],
             hi  = summ[paste0("beta_", pred_names), "hi95"]))
rownames(cmp) <- NULL
print(cmp, digits = 3)

# Station baselines: posterior of alpha_j, on log scale and back-transformed.
# exp(alpha_j) = typical PM2.5 at station j on a day with average weather.
baseline <- sapply(stations, function(s) exp(draws[, paste0("alpha_", s)]))
round(t(apply(baseline, 2, quantile, c(.025, .5, .975))), 1)

# Questions OLS can't easily answer: probability one station is dirtier
mean(draws[, paste0("alpha_", stations[1])] > draws[, paste0("alpha_", stations[2])])
# % change in PM2.5 per +1 SD of wind (with uncertainty)
round(quantile(100 * (exp(draws[, "beta_Wind"]) - 1), c(.025, .5, .975)), 1)

# ---- Plot: coefficient comparison ------------------------------------------
png("figures/coef_comparison.png", 1100, 550, res = 120)
par(mar = c(4, 7, 2, 1))
mcol <- c("OLS pooled" = "grey50", "OLS + station dummies" = "#d1603d",
          "Hierarchical (HMC)" = "#2a78b5")
ypos <- rep(1:3, each = 1) + rep(c(-0.2, 0, 0.2), each = 3)
plot(cmp$est, ypos, xlim = range(c(cmp$lo, cmp$hi)), pch = 19, col = mcol[cmp$model],
     yaxt = "n", ylab = "", xlab = "coefficient (log PM2.5 per 1 SD)",
     main = "Weather effects: OLS vs hierarchical HMC")
segments(cmp$lo, ypos, cmp$hi, ypos, col = mcol[cmp$model], lwd = 2)
axis(2, at = 1:3, labels = pred_names, las = 1); abline(v = 0, lty = 3)
legend("bottomright", names(mcol), col = mcol, pch = 19, lwd = 2, bty = "n", cex = .8)
dev.off()

# ---- Plot: station baselines ------------------------------------------------
png("figures/station_baselines.png", 1000, 550, res = 120)
scol <- c("#2a78b5", "#d1603d", "#3a9d5d")
dens <- apply(baseline, 2, density)
plot(NULL, xlim = range(sapply(dens, `[[`, "x")), ylim = c(0, max(sapply(dens, function(z) max(z$y)))),
     xlab = "Baseline PM2.5 (ug/m^3) on an average-weather day", ylab = "posterior density",
     main = "Station-specific baselines, exp(alpha_j)")
for (j in 1:J) lines(dens[[j]], col = scol[j], lwd = 3)
legend("topright", stations, col = scol, lwd = 3, bty = "n")
dev.off()

# =============================================================================
# PART 5: POSTERIOR PREDICTIVE CHECK
# =============================================================================
# Simulate new datasets from posterior draws and compare with the observed data.
png("figures/posterior_predictive.png", 1400, 480, res = 120)
par(mfrow = c(1, 3))
ids <- sample(nrow(draws), 50)
for (j in 1:J) {
  rows <- sid == j
  plot(density(y[rows]), lwd = 3, main = paste("Posterior predictive:", stations[j]),
       xlab = "log PM2.5", ylim = c(0, 0.6))
  for (k in ids) {
    mu_k <- draws[k, paste0("alpha_", stations[j])] + X[rows, ] %*% draws[k, paste0("beta_", pred_names)]
    lines(density(rnorm(sum(rows), mu_k, draws[k, "sigma"])), col = adjustcolor(scol[j], .15))
  }
  lines(density(y[rows]), lwd = 3)
}
dev.off()
