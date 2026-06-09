# =============================================================================
# 02_growth_mortality.R
# Calculate net growth rates and mortality (loss) rates.
# Requires: `data` produced by 01_preprocessing.R
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Net growth rate
# μ_net(t) = Δ ln(Ñ) / Δt   where Ñ is the smoothed abundance
# The first time point in each series receives NA (no prior observation).
# -----------------------------------------------------------------------------
data <- data %>%
  group_by(population, series, .drop = FALSE) %>%
  mutate(
    net = c(NA, diff(log(abundance_smooth)) / diff(time))
  )

# -----------------------------------------------------------------------------
# 2. Mortality (loss) rate
# Loss = μ_net(Control) − μ_net(treatment), clamped to [0, ∞).
# Calculated for bacterial and Prochlorococcus populations only;
# other populations receive NaN.
# -----------------------------------------------------------------------------
bacterial_populations <- c("Prochlorococcus", "Bacteria", "smallBacteria", "largeBacteria")

data <- data %>%
  group_by(population, replicate, time, .drop = FALSE) %>%
  mutate(
    loss = case_when(
      population %in% bacterial_populations ~ net[treatment == "Control"] - net,
      TRUE ~ NaN
    ),
    # Clamp negative loss values to zero
    loss = if_else(loss < 0, 0, loss)
  )

message("Growth and mortality rates calculated. `data` object updated.")
