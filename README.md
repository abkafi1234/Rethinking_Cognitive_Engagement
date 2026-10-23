# Critical engagement with generative AI and AI dependency: survey data and R code

This repository holds the data and analysis code for a cross-sectional survey of undergraduate computer science and engineering (CSE) students in Bangladesh.

**The question.** AI-literacy frameworks implicitly assume that critical engagement protects students against AI dependency. Critical engagement here means building on AI output and cross-checking it. The study tests that assumption: is critical engagement associated with *less* AI dependency once frequency of use and AI-use experience are taken into account?

**The answer.** No. Critical engagement was positively associated with AI dependency (standardised β = .46, 95% CI [.33, .59]).

`sem_analysis.R` reproduces every analysis reported in the article and in its supplementary material. The same script accompanies the article as Online Resource 1.

## Contents

| File | Description |
|---|---|
| `Dataset.csv` | Survey responses as exported from Google Forms: 557 respondents and 22 columns (questionnaire items 4–24 and self-reported CGPA) |
| `cleaned_constructs.csv` | The 16 items used in the models for the 552 respondents with complete data, renamed and coded numerically. The script excludes 19 straight-lining respondents, so the analysis uses 533 |
| `sem_analysis.R` | The complete analysis: data preparation, measurement model, common-method checks, hypothesis test, robustness and exploratory analyses, and all tables |
| `consent_form.docx` | The consent form shown to participants before the questionnaire |
| `ethics_approval_IEC-FSIT-DIU-2025-2008.pdf` | The ethics approval letter |
| `LICENSE` | MIT License |

## Study summary

- **Design:** cross-sectional online survey (Google Forms), 8–13 January 2026; convenience sample.
- **Participants:** undergraduate CSE students at three universities in Bangladesh. Of 557 responses, 552 were complete; the other five had no answer to item 14. Nineteen of the 552 gave identical answers to all 15 seven-point items and were excluded as careless respondents, so the analysis uses 533.
- **Measures:** seven-point scales on which higher scores mean stronger agreement or greater reliance. The items were adapted from published instruments.

  | Variable | Items | Role |
  |---|---|---|
  | Critical engagement | Q8, Q9 | Predictor |
  | AI dependency | Q6, Q7 | Outcome |
  | Acceptance without verification | Q10 | Not in the dependency scale; robustness checks R3, R9 and R10 |
  | Frequency of AI use | Q5 (1–5) | Control |
  | AI-use experience | Q4 (1–4) | Control |
  | Perceived benefits | Q17–Q21 | Covariate in a robustness check and in the exploratory analysis |
  | Self-reported CGPA | — | Covariate in a robustness check |
  | AI-induced remorse and concern | Q11, Q12, Q14–Q16 | Marker variable for the common-method checks |

- **Analysis:** covariance-based structural equation modelling in lavaan, with robust maximum likelihood. Critical engagement enters as a single indicator (the mean of Q8 and Q9) with its error variance fixed at (1 − ω) × its variance. The analysis includes:
  - bias-corrected bootstrap confidence intervals (5,000 resamples), and a two-step bootstrap that re-estimates ω in every resample;
  - ten robustness checks;
  - common-method checks, with the Lindell–Whitney marker adjustment as the primary mitigation;
  - a sensitivity analysis that fixes the reliability of critical engagement at .50 to .80;
  - an exploratory model of the indirect association through perceived benefits.

## Data

### `Dataset.csv`

Each column header is the question as worded in the questionnaire.

| Column | Item | Coding |
|---|---|---|
| 1 | Q4: AI tool usage experience | Text: `Less than 6 months`, `6-12 months`, `1-2 years`, `More than 2 years` (coded 1–4 by the script) |
| 2 | Q5: How often do you use AI tools to find information or solve problems? | Text: `Never`, `Rarely`, `Sometimes`, `Often`, `Always` (coded 1–5) |
| 3–21 | Q6–Q24: attitudinal items (the header of Q6 has no item number) | 1–7 |
| 22 | Self-reported cumulative GPA (CGPA) | Free text on a 0–4 scale; 550 values are usable |

- **Missing values:** Q14 is blank for five respondents, who are excluded from the analysis.
- **Not in the public file:** the demographic questions and the open-ended comments.
- **Privacy:** no personally identifying information was collected.

### `cleaned_constructs.csv`

The file has 552 rows. They are the complete responses of `Dataset.csv` in the same order, and the script relies on that order to add experience and CGPA. Values are stored as decimals (e.g. `4.0`). Column names combine a construct prefix, the item number and a short label:

| Prefix | Columns | Use in the analysis |
|---|---|---|
| `CE_` | Q8 (building on AI output as a starting point), Q9 (cross-checking) | Critical engagement |
| `AD_` | Q6, Q7, Q10 | AI dependency (Q6, Q7). Q10, acceptance without verification, is analysed separately (R3, R9 and R10) |
| `AD_Q5_frequency` | Q5 (1–5) | Frequency of use, a control variable. The prefix reflects the four-item dependency scale of robustness check R3, which includes this item |
| `PB_` | Q17–Q21 | Perceived benefits |
| `AR_` | Q11, Q12, Q14–Q16 | Remorse and concern (marker variable) |

The script reads AI-use experience (Q4) and CGPA from `Dataset.csv`.

## Reproducing the analysis

**Requirements.** The article used R 4.6.1 with lavaan 0.7-3, semTools 0.5-10 and psych 2.6.9. Other versions may change results in the third decimal.

```r
install.packages(c("lavaan", "semTools", "psych"))
```

**Running it.** From the folder that contains the three files, run:

```bash
Rscript sem_analysis.R
```

The run takes about 5 minutes: four bootstraps of 5,000 resamples, one of which refits the measurement model in every resample. The seed is set to 20261018. Results are written to `analysis/output/`, which is created automatically:

| File | Contents |
|---|---|
| `sem_analysis_log.txt` | The full log, Sections 0–12 |
| `manuscript_values.md` | Every number as printed in the article's tables, each rounded once from the unrounded estimate |
| `results_summary.md` | Key results to three decimals |
| `table_*.csv` | Descriptives, correlations, loadings, reliability and validity, common-method checks, the hypothesis model, the robustness checks and the reliability sensitivity |
| `online_resource_2_tables.tex` | The tables of the supplementary material (Online Resource 2) |
| `fig_S1_reliability_sensitivity.pdf` | Fig. S1 of Online Resource 2: the coefficient of critical engagement for assumed reliabilities of .50 to .80 |

Running the script prints a note from semTools, which rescales the latent variables for the discriminant-validity test, and lavaan warnings about non-admissible bootstrap draws. Both are expected. Most of the non-admissible draws come from robustness check R4, which the article uses to show why critical engagement is not modelled as a two-indicator factor.

## Structure of the script

| Section | Content | Where it appears in the article |
|---|---|---|
| 0–2 | Software versions, data preparation, descriptive statistics | Sample description and Table 3 |
| 3–4 | Measurement model, reliability and validity | Table 2 |
| 5 | Common-method variance | Common-method assessment |
| 6 | Hypothesis test, with bootstrap confidence intervals and the two-step bootstrap | Table 4 |
| 7 | Ten robustness checks | Table 5 and Fig. 2 |
| 8 | Sensitivity analyses: smallest detectable correlation; reliability assumed for critical engagement | Method; Section 4.4; Online Resource 2, Table S9 and Fig. S1 |
| 9 | Exploratory indirect association through perceived benefits | Results; Online Resource 2, Table S8 |
| 10–12 | Summary, values as printed, supplementary tables | Online Resource 2 |

## Ethics and consent

The study was approved by the Research Ethics Committee of the Faculty of Science and Information Technology, Daffodil International University (ref. IEC-FSIT/DIU/2025/2008, 7 January 2026). Participation was voluntary and anonymous, and all participants gave electronic informed consent before starting the questionnaire.

## License and citation

The repository is released under the MIT License (see `LICENSE`). If you use the data or code, please cite the article; its citation details will be added here when it is published.

## Contact

For questions about the data or code, please open an issue in this repository.
