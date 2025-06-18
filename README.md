![R-CMD-check](https://github.com/bips-hb/IDEFICS_scalc/actions/workflows/r.yml/badge.svg)
# 📦 IDEFICS.scalc

**IDEFICS.scalc** provides functions to compute standardized percentiles, z-scores, and a composite Metabolic Syndrome (MetS) score for children's anthropometric and metabolic parameters, based on the IDEFICS reference study. The package also supports categorizing health risk levels using action thresholds.

## 🔧 Installation

You can install the development version of `IDEFICS.scalc` from GitHub using:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install directly from GitHub
devtools::install_github("bips-hb/IDEFICS.scalc")
```

## 🚀 Example

```r
library(IDEFICS.scalc)

# Input: data frame with raw values
df <- data.frame(
  sex = c("f", "m"),
  age = c(6, 7),
  height = c(120, 125),
  waist = c(55, 60),
  homa = c(1.2, 1.4),
  sbp = c(100, 105),
  dbp = c(65, 70),
  trg = c(0.9, 1.0),
  hdl = c(1.1, 1.0)
)

# Calculate z-scores and MetS
results <- ScoreCalc(df, return_values = c("z.score", "MetS"))

# View output
print(results)
```

## 📚 Functions

- `get_scores()` — Get percentiles or z-scores for a single variable.
- `ScoreCalc()` — Calculate scores for multiple variables and optionally MetS.
- `MetSScore()` — Calculate MetS from z-scores.
- `action_levels()` — Assign action levels based on percentile cutoffs.

## 📖 Reference

This package uses internal parameter tables based on the **IDEFICS study**, a European cohort focused on childhood obesity and metabolic health:

-  Ahrens, W., Moreno, L. A., Mårild, S., Molnár, D., Siani, A., De Henauw, S., Böhmann, J., Günther, K., Hadjigeorgiou, C., Iacoviello, L., Lissner, L., Veidebaum, T., Pohlabeln, H. & Pigeot, I. on behalf of the IDEFICS consortium (2014). [Metabolic syndrome in young children: definitions and results of the IDEFICS study.](https://doi.org/10.1038/ijo.2014.130) *International Journal of Obesity*, 38 (Suppl 2), S4–S14.

The reference data used in this package was originally published in:
- [bips-hb/IDEFICS-Score_Calculator](https://github.com/bips-hb/IDEFICS-Score_Calculator/) — a repository containing the IDEFICS tables and an R script version of the scoring logic.

This package improves upon and packages that logic cleanly for programmatic and research use.

