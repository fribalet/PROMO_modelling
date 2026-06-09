# =============================================================================
# 03_model_fitting.R
# Fit Holling Type I and II functional response models to grazer and phage data.
# Computes goodness-of-fit metrics and bootstrap 95% CIs.
# Requires: `data` produced by 02_growth_mortality.R
# =============================================================================

library(DEoptim)

# -----------------------------------------------------------------------------
# 1. Reshape to wide format for model fitting
# Abundance is rescaled to units of 10^7 cells mL⁻¹ for numerical stability.
# -----------------------------------------------------------------------------
data_w <- data %>%
  mutate(abundance = abundance / 1e7) %>%
  pivot_wider(
    id_cols     = c(time, treatment, experiment, series),
    names_from  = population,
    values_from = c(abundance, net, loss),
    values_fill = 0
  )

# -----------------------------------------------------------------------------
# 2. Functional response models
# -----------------------------------------------------------------------------

#' Holling Type I: f(x) = a * x
model_pp <- function(a, data) data$x * a

#' Holling Type II: f(x) = (a * x) / (1 + a * h * x)
model_pp2 <- function(a, h, data) (data$x * a) / (1 + data$x * a * h)

# RMSE objective functions (reference `mydata` in calling environment)
measure_distance  <- function(params, data) sqrt(mean((data$y - model_pp( params[1],           data))^2))
measure_distance2 <- function(params, data) sqrt(mean((data$y - model_pp2(params[1], params[2], data))^2))

data_dist  <- function(params, ...) measure_distance( params, mydata)
data_dist2 <- function(params, ...) measure_distance2(params, mydata)

# -----------------------------------------------------------------------------
# 3. Prepare fitting subsets
# Only retain observations where measured loss is biologically meaningful.
# -----------------------------------------------------------------------------

# Grazer subset: low/high grazer treatments; exclude near-zero loss
mydata.Grazer <- data_w %>%
  filter(treatment %in% c("Low grazer", "High grazer")) %>%
  filter(loss_Prochlorococcus > 0.05) %>%
  mutate(
    x       = abundance_Prochlorococcus,
    predator = abundance_Grazer,
    y       = x * loss_Prochlorococcus / predator
  )

# Phage subset: low/high virus treatments; also exclude unreliably high lysis
mydata.Phage <- data_w %>%
  filter(treatment %in% c("High virus", "Low virus")) %>%
  filter(loss_Prochlorococcus > 0.05) %>%
  filter(abundance_Prochlorococcus * loss_Prochlorococcus / abundance_Phage < 0.07) %>%
  mutate(
    x        = abundance_Prochlorococcus,
    predator = abundance_Phage,
    y        = x * loss_Prochlorococcus / predator
  )

# DEoptim control settings
deoptim_ctrl <- DEoptim.control(
  reltol  = 1e-6,
  trace   = 50,
  steptol = 500,
  strategy = 2,
  itermax = 1000
)

# -----------------------------------------------------------------------------
# 4. Fit models
# -----------------------------------------------------------------------------

message("Fitting Grazer Type I ...")
mydata <- mydata.Grazer
best_g1       <- DEoptim(data_dist,  lower = 0,      upper = 1000,        control = deoptim_ctrl)
params.Grazer <- best_g1$optim$bestmem

message("Fitting Grazer Type II ...")
best_g2        <- DEoptim(data_dist2, lower = c(0,0), upper = c(1000,1000), control = deoptim_ctrl)
params.Grazer2 <- best_g2$optim$bestmem

message("Fitting Phage Type I ...")
mydata <- mydata.Phage
best_p1      <- DEoptim(data_dist,  lower = 0,      upper = 1000,        control = deoptim_ctrl)
params.Phage <- best_p1$optim$bestmem

message("Fitting Phage Type II ...")
best_p2       <- DEoptim(data_dist2, lower = c(0,0), upper = c(1000,1000), control = deoptim_ctrl)
params.Phage2 <- best_p2$optim$bestmem

# -----------------------------------------------------------------------------
# 5. Goodness-of-fit (RMSE, R², AIC)
# -----------------------------------------------------------------------------

compute_gof <- function(params, data, k, type = c("I", "II")) {
  type <- match.arg(type)
  pred <- if (type == "I") model_pp(params[1], data) else model_pp2(params[1], params[2], data)
  n    <- nrow(data)
  rss  <- sum((data$y - pred)^2)
  tibble(
    RMSE = sqrt(rss / n),
    R2   = 1 - rss / sum((data$y - mean(data$y))^2),
    AIC  = n * log(rss / n) + 2 * k,
    n    = n,
    k    = k
  )
}

gof <- bind_rows(
  compute_gof(params.Grazer,  mydata.Grazer, k = 1, type = "I")  |> mutate(dataset = "Grazer", model = "Type I"),
  compute_gof(params.Grazer2, mydata.Grazer, k = 2, type = "II") |> mutate(dataset = "Grazer", model = "Type II"),
  compute_gof(params.Phage,   mydata.Phage,  k = 1, type = "I")  |> mutate(dataset = "Phage",  model = "Type I"),
  compute_gof(params.Phage2,  mydata.Phage,  k = 2, type = "II") |> mutate(dataset = "Phage",  model = "Type II")
)

print(gof)
write_csv(gof, "outputs/gof_table.csv")

# -----------------------------------------------------------------------------
# 6. Bootstrap 95% confidence intervals
# -----------------------------------------------------------------------------

bootstrap_typeI <- function(data, n_boot = 1000, seed = 42) {
  set.seed(seed)
  replicate(n_boot, {
    d <- data[sample(nrow(data), replace = TRUE), ]
    optim(
      par    = 1,
      fn     = function(a) sqrt(mean((d$y - model_pp(a, d))^2)),
      method = "Brent", lower = 0, upper = 1000
    )$par
  })
}

bootstrap_typeII <- function(data, init_params, n_boot = 1000, seed = 42) {
  set.seed(seed)
  t(replicate(n_boot, {
    d <- data[sample(nrow(data), replace = TRUE), ]
    optim(
      par    = init_params,
      fn     = function(p) sqrt(mean((d$y - model_pp2(p[1], p[2], d))^2)),
      method = "L-BFGS-B", lower = c(0, 0), upper = c(1000, 1000)
    )$par
  }))
}

boot_grazer  <- bootstrap_typeI(mydata.Grazer)
boot_phage   <- bootstrap_typeI(mydata.Phage)
boot_grazer2 <- bootstrap_typeII(mydata.Grazer, init_params = params.Grazer2)

ci_grazer    <- quantile(boot_grazer,       c(0.025, 0.975))
ci_phage     <- quantile(boot_phage,        c(0.025, 0.975))
ci_a_grazer2 <- quantile(boot_grazer2[, 1], c(0.025, 0.975))
ci_h_grazer2 <- quantile(boot_grazer2[, 2], c(0.025, 0.975))

cat(sprintf(
  "Grazer (Type I)  a = %.4f  95%% CI [%.4f, %.4f]\n",
  params.Grazer[1], ci_grazer[1], ci_grazer[2]
))
cat(sprintf(
  "Phage  (Type I)  a = %.6f  95%% CI [%.6f, %.6f]\n",
  params.Phage[1], ci_phage[1], ci_phage[2]
))
cat(sprintf(
  "Grazer (Type II) a = %.4f  95%% CI [%.4f, %.4f]\n",
  params.Grazer2[1], ci_a_grazer2[1], ci_a_grazer2[2]
))
cat(sprintf(
  "Grazer (Type II) h = %.5f  95%% CI [%.5f, %.5f]\n",
  params.Grazer2[2], ci_h_grazer2[1], ci_h_grazer2[2]
))

# Save parameter table
param_table <- tibble(
  dataset  = c("Grazer", "Grazer", "Grazer", "Phage"),
  model    = c("Type I", "Type II", "Type II", "Type I"),
  param    = c("a", "a", "h", "a"),
  estimate = c(params.Grazer[1], params.Grazer2[1], params.Grazer2[2], params.Phage[1]),
  ci_lower = c(ci_grazer[1],    ci_a_grazer2[1],   ci_h_grazer2[1],  ci_phage[1]),
  ci_upper = c(ci_grazer[2],    ci_a_grazer2[2],   ci_h_grazer2[2],  ci_phage[2])
)
write_csv(param_table, "outputs/parameter_estimates.csv")

# -----------------------------------------------------------------------------
# 7. Theoretical mortality predictions (full dataset)
# -----------------------------------------------------------------------------
data_w3 <- data_w %>%
  ungroup() %>%
  mutate(
    x       = abundance_Prochlorococcus,
    predator = abundance_Grazer,
    y       = x * loss_Prochlorococcus / predator) %>%
  mutate(
    theoretical.loss.Phage_Prochlorococcus  = model_pp(params.Phage[1],  .) * abundance_Phage  / abundance_Prochlorococcus,
    theoretical.loss.Grazer_Prochlorococcus = model_pp(params.Grazer[1], .) * abundance_Grazer / abundance_Prochlorococcus,
    theoretical.loss_Prochlorococcus        = theoretical.loss.Grazer_Prochlorococcus +
                                               theoretical.loss.Phage_Prochlorococcus
  ) %>%
  mutate(across(contains("theoretical"), ~ if_else(. == 0, NA_real_, .))) %>%
  mutate(across(contains("abundance"),   ~ .x * 1e7)) %>%    # restore original units
  mutate(
    theoretical.loss_Prochlorococcus = if_else(
      str_detect(treatment, "virus")  & abundance_Phage  == 0, NA_real_, theoretical.loss_Prochlorococcus),
    theoretical.loss_Prochlorococcus = if_else(
      str_detect(treatment, "grazer") & abundance_Grazer == 0, NA_real_, theoretical.loss_Prochlorococcus)
  ) %>%
  select(-x, -predator, -y) %>%
  arrange(time)

message("Model fitting complete. GOF and parameter tables written to outputs/.")
