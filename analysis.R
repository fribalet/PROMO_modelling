# =============================================================================
# Master analysis script
# Prochlorococcus Mortality: Grazing and Viral Lysis
# =============================================================================
# Run this script from the repository root to execute the full analysis.
# Outputs are written to figures/ and outputs/.

# -- Setup --------------------------------------------------------------------
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("outputs")) dir.create("outputs")

# -- Run steps ----------------------------------------------------------------
source("R/01_preprocessing.R")
source("R/02_growth_mortality.R")
source("R/03_model_fitting.R")
source("R/04_figures.R")

message("Analysis complete. Check figures/ and outputs/ for results.")
