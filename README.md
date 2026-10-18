# Critical engagement and generative AI dependency: a test of the protective-factor assumption

Manuscript and analysis for a cross-sectional survey of 552 computer science and engineering undergraduates at three universities in Bangladesh. The research question: is critical engagement with generative AI associated with *less* AI dependency, as AI-literacy frameworks assume, once frequency of use and AI-use experience are controlled? It is not. Critical engagement is associated with *more* dependency (β = .46, 95% CI [.32, .60]).

- **Target journal:** Education and Information Technologies (Springer).
- **Status (2026-10-18):** submission draft. Open items for the authors are listed in [docs/revision_summary.md](docs/revision_summary.md).

## Folders

| Folder | Contents |
|---|---|
| [manuscript/](manuscript/) | `main.tex` and `main.pdf` (the paper), the anonymised version, title page, cover letter and the upload set in `submission/`; `references.bib`, figures, reference scripts; `archive/` holds the earlier Elsevier version ("Update 2") |
| [analysis/](analysis/) | `sem_analysis.R` produces every result in the paper; `output/` holds the log, tables and `manuscript_values.md` |
| [data/](data/) | `raw/Dataset.csv` (557 responses) and `processed/cleaned_constructs.csv` (552 analytic cases) |
| [docs/](docs/) | Revision summary, statistical corrections, novelty and positioning memo, literature comparison table, title and structure plan, benchmark against ten EAIT papers |
| [literature_search/](literature_search/) | All Consensus searches (S1–S14) and their results; `eait_exemplars/` holds 12 EAIT papers read for title and structure conventions |
| [journal_selection/](journal_selection/) | Journal recommendation report, ACU eligibility list and metrics |
| [ethics/](ethics/) | Ethics approval, consent form, data usage agreement. The approval for this study is ref. IEC-FSIT/DIU/2025/2008; `REC Approval.pdf` belongs to a different study (see `ethics/README.md`) |
| [github_repo/](github_repo/) | README for the public data repository (github.com/abkafi1234/Rethinking_Cognitive_Engagement) |

## Reproduce

From this folder:

```bash
Rscript analysis/sem_analysis.R          # all results, about 75 s (two sets of 5,000 bootstrap resamples)
Rscript analysis/make_figures.R          # Figure 2
cd manuscript && latexmk -pdf -outdir=build main.tex && cp build/main.pdf main.pdf
```

Requirements:

- R 4.6 with lavaan 0.7-3, semTools 0.5-10 and psych 2.6.9.
- TeX Live with apacite, newtx and TikZ.
