# Critical Engagement and AI Dependency Among University Students

A cross-sectional structural equation modeling study examining the associations
among critical engagement, perceived benefits, AI dependency, and AI-induced
remorse in undergraduate students.

## Overview

This repository contains the data, analysis code, and manuscript sections for
a study testing whether critical engagement with AI tools is associated with
*lower* AI dependency, as commonly assumed in AI literacy frameworks, or
whether that protective assumption fails to hold. Using survey data from 552
university students and covariance-based structural equation modeling
(CB-SEM), the study finds that critical engagement is **positively**
associated with AI dependency, both directly and indirectly through perceived
benefits, contrary to the protective-factor assumption found in AI literacy,
media literacy, and health behaviour research. AI-induced remorse, in turn,
is associated with perceived benefits rather than with dependency itself.

## Constructs

| Code | Construct | Items | Description |
|------|-----------|-------|--------------|
| AD | AI Dependency | 4 | Reliance on AI for decisions, problem-solving, and verification under time pressure |
| AR | AI-Induced Remorse | 5 | Guilt and regret associated with AI reliance |
| PB | Perceived Benefits | 5 | Beliefs about AI's value for time savings, confidence, and learning |
| CE | Critical Engagement | 2 | Reflective and evaluative use of AI outputs |


> **Note on file versions:** several `.R` and `.tex` files in this repository
> represent earlier iterations kept for traceability. **`full_analysis.R`** is
> the authoritative analysis script. For the manuscript text, the most recent
> file in each section (by filename and date) supersedes earlier drafts —
> see the "Which file is final" table below.

## Which File Is Final

| Section | Use this file |
|---|---|
| Data cleaning + construct scoring | `cleaned_constructs.csv` |
| Full CFA / SEM / mediation analysis | `full_analysis.R` |
| Preliminary diagnostics (KMO, Bartlett, Harman, Mardia) | `diagnostic_tests.R` (now folded into `full_analysis.R`) |
| Abstract | `abstract_theme_final.tex` |
| Introduction & Theory | `intro_theory_final.tex` |
| Methodology | `methodology_rewrite.tex` |
| Results | `results_rewrite.tex` |
| Discussion | `discussion_rewrite.tex` |
| Limitations & Conclusion | `limitations_conclusion_rewrite.tex` |

## Method Summary

- **Design:** Cross-sectional online survey (Google Forms), convenience sampling
- **Sample:** 552 undergraduate CSE students (after listwise deletion from 557 responses), three Bangladeshi universities
- **Scale:** 7-point Likert (1 = strongly disagree to 7 = strongly agree); one frequency item recoded from ordinal to numeric
- **Estimator:** MLR (Satorra–Bentler robust correction) for CFA/SEM fit indices; ML with 5,000-resample bootstrap and BCa 95% CIs for structural paths and mediation
- **Software:** R, `lavaan`, `semTools`, `psych`

## Key Results

| Path | β (std.) | p | Decision |
|------|----------|---|----------|
| AD → AR (H1) | −.032 | .749 | Not supported |
| PB → AD (H2) | .391 | <.001 | Supported |
| PB → AR (H3) | .208 | .013 | Supported |
| CE → PB (H4) | .697 | <.001 | Supported |
| CE → AD (H5) | .310 | .007 | Supported |
| CE → PB → AD, indirect (H6) | .273 | <.001 | Supported (partial mediation) |

**Model fit:** CFI = .925, TLI = .903, RMSEA = .070 [.062, .079], SRMR = .071

Full statistics, including standard errors, z-values, and 95% bias-corrected
accelerated bootstrap confidence intervals, are reported in `results_rewrite.tex`.

## Reproducing the Analysis

```r
# 1. Install dependencies
install.packages(c("lavaan", "semTools", "semPlot", "psych", "dplyr"))

# 2. Set working directory to wherever cleaned_constructs.csv is saved
setwd("path/to/this/repository")

# 3. Run the full analysis (saves output to console; redirect with sink() if needed)
source("full_analysis.R", encoding = "UTF-8")
```

Outputs include CFA fit indices, standardized loadings, reliability and
validity tables (α, CR/ω, AVE, HTMT), structural path coefficients with
bootstrapped CIs, mediation results, and a path diagram (`sem_path_diagram_final.pdf`).

## Data

`cleaned_constructs.csv` contains 552 rows and 16 item-level columns
(4 per construct, except CE with 2), recoded and renamed from the raw
`cleaned.csv` survey export. Item wording and codes are listed in the
Measures table within `methodology_rewrite.tex`.

## Limitations

- Cross-sectional design: associations, not causal effects
- Critical Engagement measured with 2 items only (tau-equivalence constraint)
- AI Dependency includes one ordinal-scale item (frequency), introducing minor scale heterogeneity
- Convenience sample from three Bangladeshi universities; generalizability untested

See `limitations_conclusion_rewrite.tex` for the full discussion.

## Citation

Citation details to be added upon acceptance/publication.

## License

MIT
