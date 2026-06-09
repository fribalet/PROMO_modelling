# =============================================================================
# 01_preprocessing.R
# Load raw data, interpolate missing values, and smooth abundance time series.
# =============================================================================

library(tidyverse)
library(zoo)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
pre_data <- read_csv("data/data.csv", show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 2. Interpolate across NAs
# NAs are filled by extending the nearest non-NA edge values within each
# population × series group.  Rows that remain NA (e.g., leading/trailing NAs
# that cannot be extended) are dropped.
# Series C2 is truncated at time < 30 due to poor extrapolation behaviour.
# -----------------------------------------------------------------------------
data <- pre_data %>%
  group_by(population, series, .drop = FALSE) %>%
  mutate(abundance = zoo::na.fill(abundance, fill = "extend")) %>%
  filter(!is.na(abundance)) %>%
  filter(!(series == "C2" & time < 30))

# -----------------------------------------------------------------------------
# 3. Smooth abundance using a cubic spline (spar = 0.5)
# A smoothed time series reduces measurement noise and stabilises derivative
# estimates used in subsequent growth-rate calculations.
# Negative smoothed values are reflected to positive (abs).
# -----------------------------------------------------------------------------
data <- data %>%
  group_by(population, series, .drop = FALSE) %>%
  mutate(
    abundance_smooth = abs(
      smooth.spline(time, abundance, spar = 0.5, all.knots = TRUE)$y
    )
  )

message("Preprocessing complete. `data` object ready.")
