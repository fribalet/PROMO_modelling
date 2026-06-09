# Data Dictionary

## File: `data.csv`

This CSV file contains time-series abundance data for microbial populations across experimental treatments.

---

### Required Columns

| Column       | Type      | Units / Values                             | Description                                                                 |
|--------------|-----------|--------------------------------------------|-----------------------------------------------------------------------------|
| `time`       | numeric   | hours (h)                                  | Time point of the measurement                                               |
| `population` | character | e.g., `"Prochlorococcus"`, `"Bacteria"`, `"smallBacteria"`, `"largeBacteria"`, `"Grazer"`, `"Phage"` | Biological population measured                      |
| `abundance`  | numeric   | cells mL⁻¹ (or equivalent)                | Measured abundance; `NA` values are allowed and will be interpolated        |
| `treatment`  | character | e.g., `"Control"`, `"Low grazer"`, `"High grazer"`, `"Low virus"`, `"High virus"`, `"Low grazer - high virus"` | Experimental treatment applied |
| `series`     | character | e.g., `"C1"`, `"C2"`, `"A4"`, …           | Unique identifier for a time series (combination of experiment + replicate) |
| `replicate`  | character | e.g., `"R1"`, `"R2"`, …                   | Replicate label within a treatment                                          |
| `experiment` | character | e.g., `"Exp1"`, `"Exp2"`, …               | Experiment identifier                                                       |

---

### Notes

- **Missing values**: `NA` entries in `abundance` are filled by edge extension (`zoo::na.fill(..., fill = "extend")`). Rows where abundance remains `NA` after interpolation are removed.
- **Known data quality issue**: Series `C2` has unreliable extrapolated values before time point 30 and is excluded (`time < 30`) during preprocessing.
- **Units**: Internally, *Prochlorococcus* and bacterial abundances are rescaled to units of $$10^7$$ cells mL⁻¹ during model fitting, then restored to original units for output.

---

### Example rows

```
time,population,abundance,treatment,series,replicate,experiment
0,Prochlorococcus,1.2e+07,Control,C1,R1,Exp1
4,Prochlorococcus,1.5e+07,Control,C1,R1,Exp1
8,Prochlorococcus,NA,Control,C1,R1,Exp1
0,Grazer,3.0e+04,High grazer,A1,R1,Exp1
...
```
