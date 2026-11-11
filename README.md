# Prochlorococcus Mortality: Grazing and Viral Lysis

This repository contains all code and data to reproduce the analysis in:

> **Virus-grazer interplay enhances virus production and particle aggregation during infection of Prochlorococcus**  
> Debbie Lindell, Michael C. G. Carlson1, Julia Weissenbach, Shay Kirzner, Sigitas Šulčius, Yotam Hulata, Gazalah Sabehi, Svetlana Goldin, Joshua S. Sacks, Karin M. Björkman, Mathilde Dugenne, Sophia Zborowsky1, Morgan D. Linney4, Stephen J. Beckett, Ran Tahan, David Demory, Angelicque E. White, Joshua S. Weitz, David M. Karl, E. Virginia Armbrust, Anitra E. Ingalls, Francois Ribalet, David A. Caron


---

## Overview

*Prochlorococcus* is the most abundant photosynthetic organism on Earth, yet the relative contributions of **grazing** (by protists) and **viral lysis** (by phages) to its mortality remain poorly constrained. This project quantifies these mortality sources using time-series dilution experiments and fits Holling Type I and Type II functional response models to grazer and phage data.

The analysis pipeline:

1. Loads and preprocesses abundance time-series data (interpolation + cubic spline smoothing)

2. Calculates per-capita **net growth rates** and **mortality (loss) rates** relative to control treatments

3. Fits **Holling Type I and II** functional response models via Differential Evolution optimization (`DEoptim`)

4. Estimates **goodness-of-fit** (RMSE, R², AIC) and **bootstrap 95% confidence intervals** on model parameters

5. Generates figures comparing observed vs. model-predicted mortality across treatments

---

## Repository Structure

```
.
├── README.md                  ← This file
├── analysis.R                 ← Master script: runs the full pipeline
├── .gitignore
├── data/
│   ├── Figure_1_SourceData.csv← Input time-series data
│   └── treatment_table.csv    ← Name of each treatment
├── R/
│   ├── 01_preprocessing.R     ← Load data, interpolate NAs, smooth with spline
│   ├── 02_growth_mortality.R  ← Net growth rates and mortality rates
│   ├── 03_model_fitting.R     ← Model fitting, GOF metrics, bootstrap CIs
│   └── 04_figures.R           ← All publication figures
├── figures/                   ← Auto-generated (not tracked by git)
└── outputs/                   ← Auto-generated CSVs (not tracked by git)
```

---

## Models

**Holling Type I** (linear, no saturation):

$$f(x) = a \cdot x$$

**Holling Type II** (saturating):

$$f(x) = \frac{a \cdot x}{1 + a \cdot h \cdot x}$$

Where $x$ is *Prochlorococcus* abundance ($10^7$ cells mL⁻¹), $a$ is the attack/clearance rate, and $h$ is the handling time.

---

## Requirements

**R** (≥ 4.1) with the following packages:

```r
install.packages(c("tidyverse", "zoo", "DEoptim"))
```

---

## How to Run

**1. Clone the repository**

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

**2. Add the data file**

Place `data.csv` in the `data/` folder. See [`data/README.md`](data/README.md) for the expected format and column definitions.

**3. Run the full analysis**

Open R (or RStudio) from the repository root and run:

```r
source("analysis.R")
```

This will execute all four steps in order and write outputs to `figures/` and `outputs/`.

**4. Run individual steps**

You can also source any step independently, as long as the preceding steps have already been run in the same R session:

```r
source("R/01_preprocessing.R")   # produces: data
source("R/02_growth_mortality.R") # produces: data (updated)
source("R/03_model_fitting.R")    # produces: model params, GOF table, bootstrap CIs
source("R/04_figures.R")          # produces: figures/*.png
```

---

## Outputs

| File | Description |
|------|-------------|
| `outputs/gof_table.csv` | RMSE, R², AIC for all model × dataset combinations |
| `outputs/parameter_estimates.csv` | Parameter estimates with 95% bootstrap CIs |
| `figures/fig_grazer_fit.png` | Holling fit for grazer data |
| `figures/fig_phage_fit.png` | Holling fit for phage data |
| `figures/fig_mortality_all.png` | Observed vs. theoretical mortality, all treatments |
| `figures/fig2_mortality_mean.png` | Mean ± SD observed vs. modelled mortality (Figure 2) |

---

## License

[MIT License](LICENSE)

## Acknowledgements

The documentation of this repository was developed with the assistance of Claude Sonnet 4.6 (Anthropic). All scientific content, model decisions, and code were authored by François Ribalet

## Contact

François Ribalet · ribalet@uw.edu · University of Washington · https://www.ribaletlab.org
