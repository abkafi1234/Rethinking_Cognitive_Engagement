# =====================================================================================
# Online Resource 1 (Supplementary Information): R script that reproduces every analysis reported in
#   "Empirical evidence challenging the protective role of critical engagement against generative AI dependency",
#   Education and Information Technologies, and in Online Resource 2.
#
# Exploratory analysis (Section 9; Online Resource 2, Table S8): the indirect association through
#   perceived benefits. It is not a hypothesis test.
#
# Research question: Is critical engagement with generative AI associated with lower AI
#   dependency after controlling for frequency of use and AI-use experience?
# H1: Critical engagement is positively associated with AI dependency, net of frequency
#     of use and AI-use experience.
#
# Primary specification:
#   - Sample: the 533 complete responses left after excluding 19 straight-lining respondents
#     (identical answers on all 15 seven-point items). R6 re-estimates the model on all 552.
#   - AI dependency: latent factor measured by Q6 and Q7. Q10 (acceptance without verification)
#     is not part of the scale: it had the lowest loading, it kept the AVE of the scale below .50,
#     and its content mirrors the verification facet of critical engagement. R9 adds it back;
#     R10 analyses it as a separate outcome.
#   - Critical engagement: single indicator (mean of Q8 and Q9) with its error variance fixed at
#     (1 - omega) x variance, omega taken from the two-factor CFA. A two-indicator latent factor
#     gives improper bootstrap solutions; it is reported as R4.
#
# Run from the project folder:   Rscript analysis/sem_analysis.R
#   or, with Dataset.csv and cleaned_constructs.csv in the same folder as this script (as in the public
#   data repository), from that folder:   Rscript sem_analysis.R
# Inputs : data/raw/Dataset.csv, data/processed/cleaned_constructs.csv (or both files in the working folder)
# Outputs: analysis/output/  (log, CSV tables, results_summary.md, manuscript_values.md,
#          online_resource_2_tables.tex, the tables of Online Resource 2, and
#          fig_S1_reliability_sensitivity.pdf, its Fig. S1)
#
# Specification note: when a latent predictor sits beside observed covariates, lavaan
# (fixed.x = TRUE) does not estimate the latent-covariate covariances, and once any of
# them is declared it also stops estimating the covariate-covariate covariance. Every
# model below therefore declares all covariances among its exogenous variables explicitly.
# =====================================================================================
suppressPackageStartupMessages({
  library(lavaan); library(semTools); library(psych)
})
set.seed(20261018)
out_dir <- "analysis/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "sem_analysis_log.txt"), split = TRUE)
r3 <- function(x) round(x, 3)
p_fmt <- function(p) ifelse(p < .001, "<.001", sprintf("%.3f", p))
d2 <- function(x) sub("^-(\\.0+)$", "\\1", sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x)))   # two decimals, no leading zero (APA), no "-.00"; round only once
section <- function(title) cat("\n\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
mk <- function(...) paste(..., sep = "\n")
exo_cov <- function(vars) {                     # all pairwise covariances among exogenous variables
  pr <- combn(vars, 2); paste(sprintf("%s ~~ %s", pr[1, ], pr[2, ]), collapse = "\n")
}

section("0. Software")
cat(R.version.string, "\nlavaan", as.character(packageVersion("lavaan")),
    "| semTools", as.character(packageVersion("semTools")),
    "| psych", as.character(packageVersion("psych")), "\n")

# -------------------------------------------------------------------------------------
section("1. Data")
in_file <- function(...) {                      # project layout first, then the working folder
  p <- c(...)[file.exists(c(...))]
  if (length(p) == 0) stop("Data file not found: ", paste(c(...), collapse = " or "), call. = FALSE)
  p[1]
}
raw <- read.csv(in_file("data/raw/Dataset.csv", "Dataset.csv"), check.names = FALSE, stringsAsFactors = FALSE)
dat <- read.csv(in_file("data/processed/cleaned_constructs.csv", "cleaned_constructs.csv"))
v14 <- raw[[grep("^14\\.", names(raw))]]
keep <- !(is.na(v14) | trimws(as.character(v14)) == "")   # the 5 excluded cases have a blank Q14
stopifnot(sum(keep) == nrow(dat))
cat("Raw responses:", nrow(raw), "| excluded (missing Q14):", sum(!keep), "| analytic N:", nrow(dat), "\n")

exp_map <- c("Less than 6 months" = 1, "6-12 months" = 2, "1-2 years" = 3, "More than 2 years" = 4)
dat$EXP  <- unname(exp_map[trimws(raw[[1]][keep])])          # Q4: AI-use experience (1-4)
dat$FREQ <- dat$AD_Q5_frequency                              # Q5: frequency of use (1-5)
cg <- suppressWarnings(as.numeric(trimws(raw[[grep("^26\\.", names(raw))]][keep])))
cg[!is.na(cg) & (cg < 0 | cg > 4)] <- NA
dat$CGPA <- cg                                              # Q26: self-reported CGPA (0-4)
stopifnot(!anyNA(dat$EXP))

ce_items <- c("CE_Q8_starting_point", "CE_Q9_crosscheck")
ad_items <- c("AD_Q6_reliance_decisions", "AD_Q7_rely_vs_independent")   # AI dependency (primary scale)
q10      <- "AD_Q10_accept_no_verify"                                    # not in the scale; used in R3, R9 and R10
pb_items <- c("PB_Q17_saves_time", "PB_Q18_confidence", "PB_Q19_learning_efficiency", "PB_Q20_satisfaction", "PB_Q21_collaborative")
ar_items <- c("AR_Q14_remorseful_avoidance", "AR_Q15_regret_dependency", "AR_Q16_regret_after_tasks", "AR_Q11_fear_weaken_thinking", "AR_Q12_aware_recall_loss")
likert15 <- c(ce_items, ad_items, q10, pb_items, ar_items)    # all seven-point items in the instrument
likert14 <- setdiff(likert15, q10)                            # the seven-point items of the four measurement models
dat$straightline <- apply(dat[, likert15], 1, function(x) length(unique(x)) == 1)
cat("Straight-lining respondents (identical answer on all 15 seven-point items):", sum(dat$straightline), "\n")
dat$CE_mean <- rowMeans(dat[, ce_items]); dat$AD_mean <- rowMeans(dat[, ad_items])
dat$PB_mean <- rowMeans(dat[, pb_items]); dat$AR_mean <- rowMeans(dat[, ar_items])
dat_all <- dat                                                # all complete responses (N = 552): robustness check R6
dat <- dat_all[!dat_all$straightline, ]                       # primary analytic sample (N = 533)
cat("Primary analytic sample (straight-lining respondents excluded): N =", nrow(dat), "\n")
dist_pct <- function(x, lab) data.frame(category = lab, n = as.vector(table(factor(x, levels = seq_along(lab)))),
                                        pct = round(100 * as.vector(table(factor(x, levels = seq_along(lab)))) / length(x), 1))
exp_dist  <- dist_pct(dat$EXP, c("less than 6 months", "6-12 months", "1-2 years", "more than 2 years"))
freq_dist <- dist_pct(dat$FREQ, c("never", "rarely", "sometimes", "often", "always"))
cat("Experience (Q4), analytic sample:\n"); print(exp_dist, row.names = FALSE)
cat("Frequency (Q5), analytic sample:\n"); print(freq_dist, row.names = FALSE)
cat("CGPA parsed for", sum(!is.na(dat$CGPA)), "of", nrow(dat), "respondents in the analytic sample\n")

# -------------------------------------------------------------------------------------
section("2. Descriptive statistics and assumptions")
d <- psych::describe(dat[, c(ce_items, ad_items, q10, "CE_mean", "AD_mean", "FREQ", "EXP", "PB_mean", "CGPA")])
desc <- data.frame(variable = rownames(d), n = d$n, mean = r3(d$mean), sd = r3(d$sd), min = d$min, max = d$max,
                   skew = r3(d$skew), kurtosis = r3(d$kurtosis))
print(desc, row.names = FALSE)
write.csv(desc, file.path(out_dir, "table_descriptives.csv"), row.names = FALSE)

cvars <- c("CE_mean", "AD_mean", "FREQ", "EXP", "PB_mean", "CGPA")
ct <- psych::corr.test(dat[, cvars], use = "pairwise")
cat("\nPearson correlations (composites; pairwise n):\n"); print(r3(ct$r))
cat("\np-values (lower triangle unadjusted):\n"); print(r3(ct$p))
write.csv(data.frame(var = rownames(ct$r), r3(ct$r)), file.path(out_dir, "table_correlations.csv"), row.names = FALSE)

mard <- psych::mardia(dat[, c(ce_items, ad_items)], plot = FALSE)
cat(sprintf("\nMardia (4 focal items): skew b1p = %.2f, chi2 = %.2f, p = %s; kurtosis b2p = %.2f, z = %.2f, p = %s\n",
            mard$b1p, mard$skew, p_fmt(mard$p.skew), mard$b2p, mard$kurtosis, p_fmt(mard$p.kurt)))
kmo <- psych::KMO(dat[, c(ce_items, ad_items)])
bart <- psych::cortest.bartlett(cor(dat[, c(ce_items, ad_items)]), n = nrow(dat))
cat(sprintf("KMO (4 focal items) = %.3f; Bartlett chi2(%d) = %.2f, p = %s\n", kmo$MSA, bart$df, bart$chisq, p_fmt(bart$p.value)))

# -------------------------------------------------------------------------------------
section("3. Measurement model (CFA, MLR)")
ad_block  <- "AD =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent"
cfa_model <- mk("CE =~ CE_Q8_starting_point + CE_Q9_crosscheck", ad_block)
fit_cfa <- cfa(cfa_model, data = dat, estimator = "MLR")
fm_names <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "chisq.scaling.factor", "cfi.robust", "tli.robust",
              "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr")
fm_cfa <- fitMeasures(fit_cfa, fm_names); print(r3(fm_cfa))
cat("Scaled test reported by lavaan for MLR:", paste(lavInspect(fit_cfa, "options")$test, collapse = " + "), "\n")
std_cfa <- standardizedSolution(fit_cfa)
load_u <- std_cfa[std_cfa$op == "=~", c("lhs", "rhs", "est.std", "se", "z", "pvalue")]   # unrounded, used below
load <- load_u; load[, 3:6] <- r3(load[, 3:6]); print(load, row.names = FALSE)
write.csv(load, file.path(out_dir, "table_loadings.csv"), row.names = FALSE)
r_ce_ad <- std_cfa[std_cfa$op == "~~" & std_cfa$lhs == "CE" & std_cfa$rhs == "AD", ]
cat(sprintf("Latent correlation CE-AD = %.3f (SE %.3f), 95%% CI [%.3f, %.3f]\n", r_ce_ad$est.std, r_ce_ad$se, r_ce_ad$ci.lower, r_ce_ad$ci.upper))

section("4. Reliability, convergent and discriminant validity")
alpha_ce <- psych::alpha(dat[, ce_items], warnings = FALSE)$total$raw_alpha
alpha_ad <- psych::alpha(dat[, ad_items], warnings = FALSE)$total$raw_alpha
r89 <- cor(dat[, ce_items])[1, 2]; r67 <- cor(dat[, ad_items])[1, 2]
sbr <- function(r) 2 * r / (1 + r)                                 # Spearman-Brown coefficient of a two-item scale
omega <- unlist(semTools::compRelSEM(fit_cfa))
lam <- split(load_u$est.std, load_u$lhs)
ave <- sapply(lam, function(l) mean(l^2))                          # Fornell & Larcker (1981): mean squared std. loading
relval <- data.frame(construct = c("CE", "AD"), items = c(2, 2), alpha = r3(c(alpha_ce, alpha_ad)),
                     spearman_brown = r3(c(sbr(r89), sbr(r67))), inter_item_r = r3(c(r89, r67)),
                     omega = r3(as.numeric(omega[c("CE", "AD")])),
                     AVE = r3(as.numeric(ave[c("CE", "AD")])), sqrt_AVE = r3(sqrt(as.numeric(ave[c("CE", "AD")]))))
print(relval, row.names = FALSE)
alpha_pb <- psych::alpha(dat[, pb_items], warnings = FALSE)$total$raw_alpha    # perceived benefits (covariate in R2)
cat(sprintf("Perceived benefits (five items): alpha = %.3f\n", alpha_pb))
cat(sprintf("Inter-item correlations: CE (Q8, Q9) = %.3f; AD (Q6, Q7) = %.3f; AVE (unrounded): CE = %.4f, AD = %.4f\n",
            r89, r67, ave["CE"], ave["AD"]))
# The same CFA with Q10 as a third dependency indicator (the reason Q10 is not part of the scale)
fit_cfa3 <- cfa(mk("CE =~ CE_Q8_starting_point + CE_Q9_crosscheck", paste(ad_block, "+", q10)), data = dat, estimator = "MLR")
st3 <- standardizedSolution(fit_cfa3); l3 <- st3[st3$op == "=~", ]
ave3 <- tapply(l3$est.std^2, l3$lhs, mean)
q10_load <- l3$est.std[l3$rhs == q10]
cat(sprintf("With Q10 as a third dependency indicator: Q10 loading = %.3f; AVE: CE = %.3f, AD = %.3f; omega AD = %.3f\n",
            q10_load, ave3["CE"], ave3["AD"], unlist(semTools::compRelSEM(fit_cfa3))["AD"]))
write.csv(relval, file.path(out_dir, "table_reliability_validity.csv"), row.names = FALSE)

htmt <- semTools::htmt(cfa_model, data = dat)
cat("\nHTMT:\n"); print(r3(htmt))
cat(sprintf("Fornell-Larcker: sqrt(AVE) CE = %.3f, AD = %.3f vs latent r = %.3f\n",
            sqrt(ave["CE"]), sqrt(ave["AD"]), r_ce_ad$est.std))
fit_free <- cfa(cfa_model, data = dat, estimator = "MLR", std.lv = TRUE)          # factor variances = 1
fit_one  <- cfa(mk(cfa_model, "CE ~~ 1*AD"), data = dat, estimator = "MLR", std.lv = TRUE)
lrt_one <- lavTestLRT(fit_free, fit_one)
cat("\nScaled chi-square difference test of r(CE, AD) = 1:\n"); print(lrt_one)
dv <- semTools::discriminantValidity(fit_cfa)
cat("\nsemTools::discriminantValidity (tests r against .90; Ronkko & Cho 2022):\n"); print(dv)

# -------------------------------------------------------------------------------------
section("5. Common-method variance diagnostics")
harman_items <- c(likert14, "AD_Q5_frequency")                 # items of the four measurement models and the frequency item
pca1 <- psych::principal(dat[, harman_items], nfactors = 1, rotate = "none")
harman_pct <- 100 * pca1$values[1] / length(harman_items)
cat(sprintf("Harman's test (%d items: the 14 seven-point items of the measurement models and the frequency item, unrotated): first component explains %.1f%% of variance\n",
            length(harman_items), harman_pct))
fit_1f <- cfa(paste("G =~", paste(c(ce_items, ad_items), collapse = " + ")), data = dat, estimator = "MLR")
cat("Single-factor CFA of the 4 focal items:\n"); print(r3(fitMeasures(fit_1f, fm_names)))

# Unmeasured latent method construct (ULMC) on the 14 seven-point items of the four measurement models
full14 <- mk("CE =~ CE_Q8_starting_point + CE_Q9_crosscheck", ad_block,
             paste("PB =~", paste(pb_items, collapse = " + ")),
             paste("AR =~", paste(ar_items, collapse = " + ")),
             "AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss")
fit_t <- cfa(full14, data = dat, estimator = "MLR")
ulmc <- mk(full14, paste("M =~", paste(paste0("m*", likert14), collapse = " + ")),
           "M ~~ 0*CE + 0*AD + 0*PB + 0*AR")
fit_m <- cfa(ulmc, data = dat, estimator = "MLR")
st_t <- standardizedSolution(fit_t); st_m <- standardizedSolution(fit_m)
r_t <- st_t[st_t$op == "~~" & st_t$lhs == "CE" & st_t$rhs == "AD", ]
r_m <- st_m[st_m$op == "~~" & st_m$lhs == "CE" & st_m$rhs == "AD", ]
m_load <- st_m[st_m$op == "=~" & st_m$lhs == "M", "est.std"]
fmt_t <- fitMeasures(fit_t, fm_names); fmt_m <- fitMeasures(fit_m, fm_names)
cat("Trait-only CFA (14 items):\n"); print(r3(fmt_t)); cat("With method factor:\n"); print(r3(fmt_m))
cat("Admissible solutions: trait-only", suppressWarnings(lavInspect(fit_t, "post.check")), "| with ULMC", suppressWarnings(lavInspect(fit_m, "post.check")), "\n")
cat(sprintf("ULMC: method loading = %.3f (share of item variance = %.3f); robust CFI change = %.3f; CE-AD r: %.3f without vs %.3f [%.3f, %.3f] with method factor\n",
            m_load[1], mean(m_load^2), fmt_m["cfi.robust"] - fmt_t["cfi.robust"], r_t$est.std, r_m$est.std, r_m$ci.lower, r_m$ci.upper))
# Marker-variable approaches. No item in the instrument is theoretically unrelated to both constructs, so
# AI-induced remorse (AR), the same-method construct least correlated with dependency, serves as the marker.
# The Lindell-Whitney adjustment is the primary statistical mitigation; the CFA marker technique is supplementary.
r_obs <- cor(dat$CE_mean, dat$AD_mean)
r_mk_ce <- cor(dat$AR_mean, dat$CE_mean); r_mk_ad <- cor(dat$AR_mean, dat$AD_mean)
r_mk  <- min(r_mk_ce, r_mk_ad)                  # conservative CMV estimate (Lindell & Whitney 2001)
r_lw  <- (r_obs - r_mk) / (1 - r_mk)
t_lw  <- r_lw / sqrt((1 - r_lw^2) / (nrow(dat) - 3)); p_lw <- 2 * pt(-abs(t_lw), df = nrow(dat) - 3)
cat(sprintf("Lindell-Whitney: observed r(CE, AD) = %.3f; marker r = %.3f; adjusted r = %.3f, t(%d) = %.2f, p = %s\n",
            r_obs, r_mk, r_lw, nrow(dat) - 3, t_lw, p_fmt(p_lw)))
cat(sprintf("Marker correlations (scale means): with CE = %.3f, with AD = %.3f; r_M is the correlation with %s\n",
            r_mk_ce, r_mk_ad, if (r_mk_ad < r_mk_ce) "AD" else "CE"))
# The same adjustment applied to every correlation among the two scale scores and the two control variables;
# the H1 regression is then re-estimated on the adjusted matrix (Malhotra, Kim & Patil 2006). Composite
# scores without correction for unreliability; normal-theory ML on the correlation matrix.
lw_vars <- c("CE_mean", "AD_mean", "FREQ", "EXP")
R_obs <- cor(dat[, lw_vars]); R_adj <- (R_obs - r_mk) / (1 - r_mk); diag(R_adj) <- 1
lw_reg <- "AD_mean ~ CE_mean + FREQ + EXP"
fit_lw0 <- sem(lw_reg, sample.cov = R_obs, sample.nobs = nrow(dat))
fit_lw1 <- sem(lw_reg, sample.cov = R_adj, sample.nobs = nrow(dat))
lw_coef <- function(fit, x = "CE_mean") { ss <- standardizedSolution(fit); ss[ss$op == "~" & ss$rhs == x, c("est.std", "ci.lower", "ci.upper", "pvalue")] }
lw_b0 <- lw_coef(fit_lw0); lw_b1 <- lw_coef(fit_lw1)
lw_r2 <- c(lavInspect(fit_lw0, "rsquare"), lavInspect(fit_lw1, "rsquare"))
cat("Lindell-Whitney-adjusted correlation matrix (scale scores and controls):\n"); print(r3(R_adj))
# A more conservative estimate of method variance: the larger of the two marker correlations
r_mk2 <- max(r_mk_ce, r_mk_ad)
r_lw2 <- (r_obs - r_mk2) / (1 - r_mk2)
t_lw2 <- r_lw2 / sqrt((1 - r_lw2^2) / (nrow(dat) - 3)); p_lw2 <- 2 * pt(-abs(t_lw2), df = nrow(dat) - 3)
R_adj2 <- (R_obs - r_mk2) / (1 - r_mk2); diag(R_adj2) <- 1
stopifnot(min(eigen(R_adj2)$values) > 0)
fit_lw2 <- sem(lw_reg, sample.cov = R_adj2, sample.nobs = nrow(dat))
lw_b2 <- lw_coef(fit_lw2)
cat(sprintf("Regression of AD on CE, FREQ and EXP (scale scores): beta(CE) = %.3f [%.3f, %.3f] unadjusted; %.3f [%.3f, %.3f], p = %s after the adjustment (%.0f%% smaller); R2 %.3f -> %.3f\n",
            lw_b0$est.std, lw_b0$ci.lower, lw_b0$ci.upper, lw_b1$est.std, lw_b1$ci.lower, lw_b1$ci.upper, p_fmt(lw_b1$pvalue),
            100 * (1 - lw_b1$est.std / lw_b0$est.std), lw_r2[1], lw_r2[2]))
cat(sprintf("Conservative variant (r_M = %.3f, the larger marker correlation): adjusted r = %.3f, t(%d) = %.2f, p = %s; scale-score beta(CE) = %.3f [%.3f, %.3f], p = %s (%.0f%% smaller)\n",
            r_mk2, r_lw2, nrow(dat) - 3, t_lw2, p_fmt(p_lw2), lw_b2$est.std, lw_b2$ci.lower, lw_b2$ci.upper, p_fmt(lw_b2$pvalue),
            100 * (1 - lw_b2$est.std / lw_b0$est.std)))
# CFA marker technique (Williams, Hartman & Cavazotte 2010): marker factor identified by its own items,
# with equal method loadings on the focal items and no correlation with the focal factors
cfa_mk <- mk(cfa_model, paste("MK =~", paste(ar_items, collapse = " + ")), "AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss",
             paste("MK =~", paste(paste0("mm*", c(ce_items, ad_items)), collapse = " + ")), "MK ~~ 0*CE + 0*AD")
fit_mk <- cfa(cfa_mk, data = dat, estimator = "MLR")
st_mk <- standardizedSolution(fit_mk)
r_cmk <- st_mk[st_mk$op == "~~" & st_mk$lhs == "CE" & st_mk$rhs == "AD", ]
mm_load <- st_mk[st_mk$op == "=~" & st_mk$lhs == "MK" & st_mk$rhs %in% c(ce_items, ad_items), "est.std"]
cat(sprintf("CFA marker technique: method loading on focal items = %.3f (%.1f%% of variance); CE-AD r = %.3f [%.3f, %.3f]; admissible: %s\n",
            mean(mm_load), 100 * mean(mm_load^2), r_cmk$est.std, r_cmk$ci.lower, r_cmk$ci.upper, suppressWarnings(lavInspect(fit_mk, "post.check"))))
cmb <- data.frame(check = c("Harman first unrotated component (% variance, 15 items)",
                            "Single-factor CFA robust CFI (4 focal items)", "Two-factor CFA robust CFI (4 focal items)",
                            "ULMC standardized method loading (14 items)", "ULMC share of item variance",
                            "Robust CFI change when method factor added",
                            "CE-AD latent r, trait-only 14-item CFA", "CE-AD latent r, with ULMC method factor",
                            "CE-AD latent r, ULMC 95% CI lower", "CE-AD latent r, ULMC 95% CI upper",
                            "Lindell-Whitney marker r (AR)", "Lindell-Whitney adjusted CE-AD r", "Lindell-Whitney p",
                            "Scale-score regression beta(CE), unadjusted", "Scale-score regression beta(CE), Lindell-Whitney adjusted",
                            "Adjusted beta(CE) 95% CI lower", "Adjusted beta(CE) 95% CI upper", "Adjusted beta(CE) p",
                            "Conservative r_M (larger marker correlation)", "Conservative adjusted CE-AD r",
                            "Conservative adjusted beta(CE)", "Conservative adjusted beta(CE) p",
                            "CFA marker technique: CE-AD latent r", "CFA marker: 95% CI lower", "CFA marker: 95% CI upper"),
                  value = c(r3(c(harman_pct, fitMeasures(fit_1f, "cfi.robust"), fm_cfa["cfi.robust"], m_load[1], mean(m_load^2),
                               fmt_m["cfi.robust"] - fmt_t["cfi.robust"], r_t$est.std, r_m$est.std, r_m$ci.lower, r_m$ci.upper,
                               r_mk, r_lw)), p_fmt(p_lw), r3(c(lw_b0$est.std, lw_b1$est.std, lw_b1$ci.lower, lw_b1$ci.upper)),
                            p_fmt(lw_b1$pvalue), r3(c(r_mk2, r_lw2, lw_b2$est.std)), p_fmt(lw_b2$pvalue),
                            r3(c(r_cmk$est.std, r_cmk$ci.lower, r_cmk$ci.upper))))
write.csv(cmb, file.path(out_dir, "table_cmb.csv"), row.names = FALSE)

# -------------------------------------------------------------------------------------
section("6. Hypothesis test (H1): AD ~ CE + FREQ + EXP")
# Critical engagement enters as a single indicator: the mean of Q8 and Q9, with its error variance fixed at
# (1 - omega) x its variance (omega from the two-factor CFA of Section 3). The latent variable CE is therefore
# corrected for unreliability without estimating a two-indicator factor, which gives improper solutions (R4).
si_ce <- function(data, rel = as.numeric(omega["CE"])) {
  mk("CE =~ 1*CE_mean", sprintf("CE_mean ~~ %.6f*CE_mean", (1 - rel) * var(data$CE_mean)))
}
ev_ce <- (1 - as.numeric(omega["CE"])) * var(dat$CE_mean)
cat(sprintf("Single indicator: omega(CE) = %.3f; variance of the CE mean = %.3f; fixed error variance = %.3f\n",
            omega["CE"], var(dat$CE_mean), ev_ce))
h1_model <- mk(si_ce(dat), ad_block, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP")))
fit_h1 <- sem(h1_model, data = dat, estimator = "MLR")
fm_h1 <- fitMeasures(fit_h1, fm_names); print(r3(fm_h1))
cat("Top modification indices:\n"); print(modindices(fit_h1, sort. = TRUE, maximum.number = 5)[, c("lhs", "op", "rhs", "mi", "sepc.all")])
pe_h1 <- parameterEstimates(fit_h1); ss_h1 <- standardizedSolution(fit_h1)
paths <- merge(pe_h1[pe_h1$op == "~", c("lhs", "op", "rhs", "est", "se", "z", "pvalue")],
               ss_h1[ss_h1$op == "~", c("lhs", "op", "rhs", "est.std", "ci.lower", "ci.upper")], by = c("lhs", "op", "rhs"))
covs <- ss_h1[ss_h1$op == "~~" & ss_h1$lhs != ss_h1$rhs & ss_h1$lhs %in% c("CE", "FREQ", "EXP"), c("lhs", "rhs", "est.std", "pvalue")]
cat("Correlations among predictors:\n"); print(r3(covs[, 3:4])); print(covs[, 1:2])
r2_h1 <- lavInspect(fit_h1, "rsquare")["AD"]
fit_h1_0 <- sem(mk(si_ce(dat), ad_block, "AD ~ 0*CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
r2_h1_0 <- lavInspect(fit_h1_0, "rsquare")["AD"]
f2_ce <- (r2_h1 - r2_h1_0) / (1 - r2_h1)
cat(sprintf("R2(AD) = %.3f; with the CE path fixed to 0, R2 = %.3f; f2 for CE = %.3f\n", r2_h1, r2_h1_0, f2_ce))
lrt_ce <- lavTestLRT(fit_h1, fit_h1_0)
cat(sprintf("Scaled chi-square difference test for the CE path: chi2(1) = %.2f, p = %s\n", lrt_ce[2, "Chisq diff"], p_fmt(lrt_ce[2, "Pr(>Chisq)"])))

cat("\nBootstrapping (5,000 resamples, ML) ...\n")
fit_h1_boot <- sem(h1_model, data = dat, estimator = "ML", se = "bootstrap", bootstrap = 5000)
pe_boot <- parameterEstimates(fit_h1_boot, boot.ci.type = "bca.simple")   # bias-corrected (BC), no acceleration
bb <- pe_boot[pe_boot$op == "~", c("lhs", "op", "rhs", "se", "ci.lower", "ci.upper")]
names(bb)[4:6] <- c("boot_se", "bc_ci_lower", "bc_ci_upper")
# lavaan drops failed draws but keeps non-admissible ones (e.g. a negative residual variance) in the CI.
# Report how many there were, which variance went negative, and the BC CI without them.
B_h1 <- lavInspect(fit_h1_boot, "boot")
fr <- parTable(fit_h1_boot); fr <- fr[fr$free > 0, ]; fr <- fr[order(fr$free), ]
colnames(B_h1) <- paste0(fr$lhs, fr$op, fr$rhs)
nonadm <- attr(B_h1, "nonadmissible"); n_nonadm <- length(nonadm)
neg_var <- colSums(B_h1[, fr$op == "~~" & fr$lhs == fr$rhs, drop = FALSE] < 0)
cat("Failed bootstrap draws:", length(attr(B_h1, "error.idx")), "of 5000 | non-admissible draws (retained, as lavaan does):", n_nonadm, "\n")
cat("Negative variance estimates among the draws:", paste(names(neg_var)[neg_var > 0], neg_var[neg_var > 0], collapse = "; "), "\n")
bc_ci <- function(x, est) {                      # same rule as lavaan's "bca.simple"
  z0 <- qnorm(mean(x < est)); quantile(x, pnorm(2 * z0 + qnorm(c(.025, .975))), names = FALSE)
}
ok <- setdiff(seq_len(nrow(B_h1)), nonadm)
bc_adm <- t(sapply(paste0(bb$lhs, "~", bb$rhs), function(p) bc_ci(B_h1[ok, p], coef(fit_h1_boot)[p])))
bb$bc_ci_lower_admissible <- bc_adm[, 1]; bb$bc_ci_upper_admissible <- bc_adm[, 2]
h1_tab <- merge(paths, bb, by = c("lhs", "op", "rhs"))
h1_apa <- with(h1_tab, data.frame(rhs, B = sprintf("%.2f", est), SE = sprintf("%.2f", se), z = sprintf("%.2f", z), p = p_fmt(pvalue),
  beta = sprintf("%s [%s, %s]", d2(est.std), d2(ci.lower), d2(ci.upper)), bc_ci = sprintf("[%.2f, %.2f]", bc_ci_lower, bc_ci_upper),
  bc_ci_admissible = sprintf("[%.2f, %.2f]", bc_ci_lower_admissible, bc_ci_upper_admissible)))
h1_tab$pvalue <- p_fmt(h1_tab$pvalue)
num <- sapply(h1_tab, is.numeric); h1_tab[num] <- lapply(h1_tab[num], r3)
print(h1_tab, row.names = FALSE)
write.csv(h1_tab, file.path(out_dir, "table_h1_structural.csv"), row.names = FALSE)

# The bootstrap above treats omega as known. Fixing an estimated error variance as if it were known makes
# standard errors too small (Oberski & Satorra 2013), so the whole two-step procedure is also bootstrapped:
# in every resample the CFA is re-estimated, omega recomputed and the primary model refitted.
cat("\nTwo-step bootstrap (5,000 resamples): omega re-estimated in every resample ...\n")
boot2 <- t(sapply(seq_len(5000), function(b) {
  db <- dat[sample(nrow(dat), replace = TRUE), ]
  fc <- tryCatch(suppressWarnings(cfa(cfa_model, data = db, estimator = "ML")), error = function(e) NULL)
  if (is.null(fc)) return(c(omega = NA, proper_cfa = NA, B = NA))
  om_b <- suppressWarnings(unlist(semTools::compRelSEM(fc))[["CE"]])
  ok_cfa <- suppressWarnings(lavInspect(fc, "post.check"))
  if (!is.finite(om_b) || om_b <= 0 || om_b >= 1) return(c(omega = om_b, proper_cfa = ok_cfa, B = NA))
  fs <- tryCatch(suppressWarnings(sem(mk(si_ce(db, om_b), ad_block, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))),
                                      data = db, estimator = "ML")), error = function(e) NULL)
  if (is.null(fs)) return(c(omega = om_b, proper_cfa = ok_cfa, B = NA))
  c(omega = om_b, proper_cfa = ok_cfa, B = unname(coef(fs)["AD~CE"]))
}))
b2_ok <- !is.na(boot2[, "B"])
bc2 <- bc_ci(boot2[b2_ok, "B"], coef(fit_h1)[["AD~CE"]])          # centred on the unrounded estimate
cat(sprintf("Two-step bootstrap: usable resamples = %d of 5000 (omega outside (0, 1) or failed fit: %d); improper CFA in %d resamples; omega 2.5%%-97.5%% = [%.3f, %.3f]; BC CI for B = [%.3f, %.3f]\n",
            sum(b2_ok), sum(!b2_ok), sum(boot2[, "proper_cfa"] == 0, na.rm = TRUE),
            quantile(boot2[, "omega"], .025, na.rm = TRUE), quantile(boot2[, "omega"], .975, na.rm = TRUE), bc2[1], bc2[2]))

# -------------------------------------------------------------------------------------
section("7. Robustness checks (CE -> AD under alternative specifications)")
get_ce <- function(fit, label, n_used, pred = "CE", dv = "AD") {
  ss <- standardizedSolution(fit); pe <- parameterEstimates(fit)
  s <- ss[ss$op == "~" & ss$lhs == dv & ss$rhs == pred, ]; e <- pe[pe$op == "~" & pe$lhs == dv & pe$rhs == pred, ]
  est <- lavInspect(fit, "options")$estimator
  fmx <- if (est == "DWLS") fitMeasures(fit, c("cfi.scaled", "rmsea.scaled", "srmr")) else fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))
  r2 <- lavInspect(fit, "rsquare"); r2 <- if (dv %in% names(r2)) r2[dv] else NA
  data.frame(check = label, predictor = pred, n = n_used, B = r3(e$est), SE = r3(e$se), p = p_fmt(e$pvalue),
             beta = r3(s$est.std), beta_ci_lower = r3(s$ci.lower), beta_ci_upper = r3(s$ci.upper),
             R2 = r3(r2), CFI = r3(fmx[1]), RMSEA = r3(fmx[2]), SRMR = r3(fmx[3]),
             admissible = suppressWarnings(lavInspect(fit, "post.check")),   # FALSE = improper solution
             beta_apa = sprintf("%s [%s, %s]", d2(s$est.std), d2(s$ci.lower), d2(s$ci.upper)), R2_apa = d2(r2))
}
rob <- list()
rob[["main"]] <- get_ce(fit_h1, "Primary model (H1): controls frequency and experience", nrow(dat))
fit_r1 <- sem(mk(si_ce(dat), ad_block, "AD ~ CE"), data = dat, estimator = "MLR")
rob[["r1"]] <- get_ce(fit_r1, "R1 No control variables", nrow(dat))
pb_block <- paste("PB =~", paste(pb_items, collapse = " + "))
fit_r2 <- sem(mk(si_ce(dat), ad_block, pb_block, "AD ~ CE + PB + FREQ + EXP", exo_cov(c("CE", "PB", "FREQ", "EXP"))), data = dat, estimator = "MLR")
rob[["r2"]] <- get_ce(fit_r2, "R2 Adds perceived benefits as a covariate", nrow(dat))
rob[["r2pb"]] <- get_ce(fit_r2, "R2 (perceived-benefits path in the same model)", nrow(dat), pred = "PB")
fit_r3 <- sem(mk(si_ce(dat), "AD4 =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent + AD_Q10_accept_no_verify + AD_Q5_frequency",
                 "AD4 ~ CE + EXP", "CE ~~ EXP"), data = dat, estimator = "MLR")
rob[["r3"]] <- get_ce(fit_r3, "R3 Original 4-item AD incl. frequency item (no frequency control)", nrow(dat), dv = "AD4")
# R4: alternative treatments of measurement error in critical engagement
#   (a) a latent factor with two indicators; its bootstrap is run to document the improper solutions
#   (b) both constructs as reliability-corrected single indicators (saturated structural model)
#   (c) the single indicator corrected with the Spearman-Brown coefficient instead of omega
h1_latent <- mk(cfa_model, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP")))
fit_r4a <- sem(h1_latent, data = dat, estimator = "MLR")
rob[["r4a"]] <- get_ce(fit_r4a, "R4 Critical engagement as a latent factor with two indicators", nrow(dat))
cat("\nBootstrapping R4a (5,000 resamples, ML) to count improper solutions ...\n")
fit_r4a_boot <- sem(h1_latent, data = dat, estimator = "ML", se = "bootstrap", bootstrap = 5000)
B_r4 <- lavInspect(fit_r4a_boot, "boot")
fr4 <- parTable(fit_r4a_boot); fr4 <- fr4[fr4$free > 0, ]; fr4 <- fr4[order(fr4$free), ]
colnames(B_r4) <- paste0(fr4$lhs, fr4$op, fr4$rhs)
n_nonadm_r4 <- length(attr(B_r4, "nonadmissible"))
neg_var_r4 <- colSums(B_r4[, fr4$op == "~~" & fr4$lhs == fr4$rhs, drop = FALSE] < 0)
cat(sprintf("R4a bootstrap: %d of 5000 draws non-admissible (%.1f%%); negative variances: %s\n", n_nonadm_r4, 100 * n_nonadm_r4 / 5000,
            paste(names(neg_var_r4)[neg_var_r4 > 0], neg_var_r4[neg_var_r4 > 0], collapse = "; ")))
si_ad <- function(data, rel = as.numeric(omega["AD"])) {
  mk("ADs =~ 1*AD_mean", sprintf("AD_mean ~~ %.6f*AD_mean", (1 - rel) * var(data$AD_mean)))
}
fit_r4b <- sem(mk(si_ce(dat), si_ad(dat), "ADs ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
rob[["r4b"]] <- get_ce(fit_r4b, "R4 Both constructs as reliability-corrected single indicators", nrow(dat), dv = "ADs")
fit_r4c <- sem(mk(si_ce(dat, sbr(r89)), ad_block, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
rob[["r4c"]] <- get_ce(fit_r4c, "R4 Single indicator corrected with the Spearman-Brown coefficient", nrow(dat))
fit_r5 <- sem(mk(ad_block, "AD ~ CE_Q8_starting_point + CE_Q9_crosscheck + FREQ + EXP"), data = dat, estimator = "MLR")
rob[["r5a"]] <- get_ce(fit_r5, "R5 Item level: 'starting point' item (Q8)", nrow(dat), pred = "CE_Q8_starting_point")
rob[["r5b"]] <- get_ce(fit_r5, "R5 Item level: 'cross-checking' item (Q9)", nrow(dat), pred = "CE_Q9_crosscheck")
fit_r5c <- sem(mk(ad_block, "AD ~ CE_Q9_crosscheck + FREQ + EXP"), data = dat, estimator = "MLR")
rob[["r5c"]] <- get_ce(fit_r5c, "R5 Item level: cross-checking item (Q9) alone", nrow(dat), pred = "CE_Q9_crosscheck")
# R6: the primary specification on all 552 complete responses (omega re-estimated in that sample)
fit_cfa_all <- cfa(cfa_model, data = dat_all, estimator = "MLR")
omega_all <- unlist(semTools::compRelSEM(fit_cfa_all))
fit_r6 <- sem(mk(si_ce(dat_all, omega_all[["CE"]]), ad_block, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))),
              data = dat_all, estimator = "MLR")
rob[["r6"]] <- get_ce(fit_r6, "R6 All respondents, including straight-lining respondents", nrow(dat_all))
cat(sprintf("R6: omega(CE) in the full sample = %.3f\n", omega_all[["CE"]]))
# R7: ordinal items require item-level indicators, so critical engagement is a latent factor here
fit_r7 <- sem(h1_latent, data = dat, estimator = "WLSMV", ordered = c(ce_items, ad_items))
rob[["r7"]] <- get_ce(fit_r7, "R7 Ordinal estimation (WLSMV, items ordered; latent critical engagement)", nrow(dat))
d_cg <- dat[!is.na(dat$CGPA), ]
fit_r8 <- sem(mk(si_ce(d_cg), ad_block, "AD ~ CE + FREQ + EXP + CGPA", exo_cov(c("CE", "FREQ", "EXP", "CGPA"))), data = d_cg, estimator = "MLR")
rob[["r8"]] <- get_ce(fit_r8, "R8 Adds self-reported CGPA as a covariate", nrow(d_cg))
fit_r9 <- sem(mk(si_ce(dat), paste("AD3 =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent +", q10),
                 "AD3 ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
rob[["r9"]] <- get_ce(fit_r9, "R9 Dependency including acceptance without verification (Q6, Q7, Q10)", nrow(dat), dv = "AD3")
fit_r10 <- sem(mk(si_ce(dat), "AD_Q10_accept_no_verify ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
rob[["r10"]] <- get_ce(fit_r10, "R10 Outcome = accepting AI output without verification (Q10 only)", nrow(dat), dv = "AD_Q10_accept_no_verify")
robust <- do.call(rbind, rob)
stopifnot(all(robust$admissible))                # every reported model must be a proper solution
print(robust, row.names = FALSE)
write.csv(robust, file.path(out_dir, "table_robustness.csv"), row.names = FALSE)
s2 <- standardizedSolution(fit_r2)
cat("\nR2 model: latent r(CE, PB) =", r3(s2[s2$op == "~~" & ((s2$lhs == "CE" & s2$rhs == "PB") | (s2$lhs == "PB" & s2$rhs == "CE")), "est.std"]), "\n")

section("8. Sensitivity analyses: detectable effect size and the reliability assumed for critical engagement")
r_min <- tanh((qnorm(.975) + qnorm(.80)) / sqrt(nrow(dat) - 3))
cat(sprintf("With N = %d, a two-tailed test at alpha = .05 has 80%% power to detect r >= %.3f\n", nrow(dat), r_min))

# The primary model fixes the error variance of critical engagement at (1 - omega) x variance. Here the primary model
# is re-estimated with the reliability fixed at .50 to .80 in steps of .05 (Online Resource 2, Table S9 and Fig. S1).
# The range covers the 95% range of omega in the two-step bootstrap of Section 6.
rel_grid <- seq(.50, .80, by = .05)
sens_rel <- do.call(rbind, lapply(rel_grid, function(rel) {
  f <- sem(mk(si_ce(dat, rel), ad_block, "AD ~ CE + FREQ + EXP", exo_cov(c("CE", "FREQ", "EXP"))), data = dat, estimator = "MLR")
  pe <- parameterEstimates(f); ss <- standardizedSolution(f)
  e <- pe[pe$op == "~" & pe$lhs == "AD" & pe$rhs == "CE", ]; s <- ss[ss$op == "~" & ss$lhs == "AD" & ss$rhs == "CE", ]
  data.frame(reliability = rel, error_variance = (1 - rel) * var(dat$CE_mean), B = e$est, SE = e$se, z = e$z, p = e$pvalue,
             beta = s$est.std, beta_ci_lower = s$ci.lower, beta_ci_upper = s$ci.upper, R2 = lavInspect(f, "rsquare")[["AD"]],
             admissible = suppressWarnings(lavInspect(f, "post.check")))
}))
stopifnot(all(sens_rel$admissible))              # every model in the grid must be a proper solution
cat("Primary model re-estimated for assumed reliabilities of critical engagement (Online Resource 2, Table S9):\n")
print(data.frame(sens_rel[1], lapply(sens_rel[2:10], r3), sens_rel[11]), row.names = FALSE)
write.csv(data.frame(lapply(sens_rel, function(x) if (is.numeric(x)) r3(x) else x)),
          file.path(out_dir, "table_reliability_sensitivity.csv"), row.names = FALSE)
cat(sprintf("Across the grid: beta %.3f to %.3f; B %.3f to %.3f; largest p = %.2e\n",
            min(sens_rel$beta), max(sens_rel$beta), min(sens_rel$B), max(sens_rel$B), max(sens_rel$p)))
# Fig. S1 of Online Resource 2: B (+/- 1.96 robust SE) and beta (delta-method 95% CI) across the grid
om_rng <- quantile(boot2[, "omega"], c(.025, .975), na.rm = TRUE)
fig_s1 <- file.path(out_dir, "fig_S1_reliability_sensitivity.pdf")
if (capabilities("cairo")) cairo_pdf(fig_s1, width = 7.2, height = 3.3, family = "Helvetica") else pdf(fig_s1, width = 7.2, height = 3.3, family = "Helvetica")
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2, 0.8), cex = 0.78)
panel <- function(est, lo, hi, ylab, title) {
  plot(rel_grid, est, type = "n", xlim = c(.48, .82), ylim = c(0, max(hi) * 1.08), xaxt = "n", bty = "n",
       xlab = "Reliability assumed for critical engagement", ylab = ylab, main = title, font.main = 1, cex.main = 1)
  rect(om_rng[1], 0, om_rng[2], max(hi) * 1.08, col = "grey92", border = NA)
  abline(h = 0, lty = 2, col = "grey40"); abline(v = omega[["CE"]], lty = 3, col = "grey30")
  segments(rel_grid, lo, rel_grid, hi, lwd = 1.6); points(rel_grid, est, pch = 16)
  axis(1, at = rel_grid, labels = sprintf("%.2f", rel_grid))
}
panel(sens_rel$B, sens_rel$B - qnorm(.975) * sens_rel$SE, sens_rel$B + qnorm(.975) * sens_rel$SE,
      "Unstandardised coefficient B (95% CI)", "(a) B with robust standard errors")
panel(sens_rel$beta, sens_rel$beta_ci_lower, sens_rel$beta_ci_upper,
      "Standardised coefficient (95% CI)", "(b) Standardised coefficient")
invisible(dev.off())
cat("Fig. S1 written to", fig_s1, "\n")

# -------------------------------------------------------------------------------------
section("9. Exploratory analysis: indirect association through perceived benefits (Online Resource 2, Table S8)")
# Not a hypothesis test. R2 is re-specified with perceived benefits (PB) as an intervening variable:
# CE -> PB (a), PB -> AD (b) and CE -> AD (c'), with frequency and experience controlled in both equations.
# The model re-parameterises R2 (identical fit and AD equation). With cross-sectional data, the indirect
# association a*b shows how much of the CE-AD association is shared with PB but cannot establish mediation:
# the model with the reverse ordering (PB -> CE -> AD) fits identically and is reported alongside it.
med_model <- mk(si_ce(dat), ad_block, pb_block, "PB ~ a*CE + FREQ + EXP", "AD ~ cp*CE + b*PB + FREQ + EXP",
                exo_cov(c("CE", "FREQ", "EXP")), "ind := a*b", "tot := cp + a*b", "prop := a*b / (cp + a*b)")
fit_med <- sem(med_model, data = dat, estimator = "MLR")
rev_model <- mk(si_ce(dat), ad_block, pb_block, "CE ~ a2*PB + FREQ + EXP", "AD ~ cp2*PB + b2*CE + FREQ + EXP",
                exo_cov(c("PB", "FREQ", "EXP")), "ind2 := a2*b2", "tot2 := cp2 + a2*b2")
fit_rev <- sem(rev_model, data = dat, estimator = "MLR")
chi_eq <- sapply(list(R2 = fit_r2, intervening = fit_med, reverse = fit_rev), function(f) unname(fitMeasures(f, "chisq")))
cat("Chi-square of R2 and of its two re-parameterisations:", paste(names(chi_eq), sprintf("%.4f", chi_eq), collapse = "; "), "\n")
stopifnot(diff(range(chi_eq)) < 1e-4, suppressWarnings(lavInspect(fit_med, "post.check")), suppressWarnings(lavInspect(fit_rev, "post.check")))
med_est <- function(fit, labels) {
  pe <- parameterEstimates(fit); ss <- standardizedSolution(fit)
  k <- match(labels, pe$label); m <- match(paste(pe$lhs, pe$op, pe$rhs)[k], paste(ss$lhs, ss$op, ss$rhs))
  data.frame(param = labels, est = pe$est[k], se = pe$se[k], z = pe$z[k], pvalue = pe$pvalue[k],
             est.std = ss$est.std[m], ci.lower = ss$ci.lower[m], ci.upper = ss$ci.upper[m])
}
me <- med_est(fit_med, c("a", "b", "cp", "ind", "tot", "prop"))
re <- med_est(fit_rev, c("a2", "b2", "cp2", "ind2", "tot2"))
r2_med <- lavInspect(fit_med, "rsquare")
cat("\nBootstrapping the exploratory model (5,000 resamples, ML) ...\n")
fit_med_boot <- sem(med_model, data = dat, estimator = "ML", se = "bootstrap", bootstrap = 5000)
pe_mb <- parameterEstimates(fit_med_boot, boot.ci.type = "bca.simple")
me$bc_lower <- pe_mb$ci.lower[match(me$param, pe_mb$label)]; me$bc_upper <- pe_mb$ci.upper[match(me$param, pe_mb$label)]
B_med <- lavInspect(fit_med_boot, "boot")
fr_m <- parTable(fit_med_boot); fr_m <- fr_m[fr_m$free > 0, ]; fr_m <- fr_m[order(fr_m$free), ]
colnames(B_med) <- ifelse(fr_m$label != "", fr_m$label, paste0(fr_m$lhs, fr_m$op, fr_m$rhs))
nonadm_med <- attr(B_med, "nonadmissible"); ok_m <- setdiff(seq_len(nrow(B_med)), nonadm_med)
ind_bc_adm <- bc_ci((B_med[, "a"] * B_med[, "b"])[ok_m], me$est[me$param == "ind"])
cat("Failed draws:", length(attr(B_med, "error.idx")), "| non-admissible draws (retained, as lavaan does):", length(nonadm_med),
    sprintf("| BC CI for the indirect association without them: [%.3f, %.3f]\n", ind_bc_adm[1], ind_bc_adm[2]))
cat("\nIntervening-variable model (a = CE -> PB, b = PB -> AD, cp = direct CE -> AD, prop = ind / tot):\n")
print(data.frame(me[1], lapply(me[-1], r3)), row.names = FALSE)
cat(sprintf("R2: PB = %.3f; AD = %.3f\n", r2_med["PB"], r2_med["AD"]))
cat("\nEquivalent model with the reverse ordering (a2 = PB -> CE, b2 = CE -> AD, cp2 = direct PB -> AD):\n")
print(data.frame(re[1], lapply(re[-1], r3)), row.names = FALSE)

# -------------------------------------------------------------------------------------
section("10. Summary for the manuscript")
h1_ce <- h1_tab[h1_tab$rhs == "CE", ]; h1_fr <- h1_tab[h1_tab$rhs == "FREQ", ]; h1_ex <- h1_tab[h1_tab$rhs == "EXP", ]
summ <- c(
  "# Key results (generated by analysis/sem_analysis.R)", "",
  sprintf("- Software: %s; lavaan %s; semTools %s; psych %s", R.version.string, packageVersion("lavaan"), packageVersion("semTools"), packageVersion("psych")),
  sprintf("- N = %d analysed: 557 collected; 5 excluded for a missing Q14; %d straight-lining respondents excluded from the %d complete responses",
          nrow(dat), sum(dat_all$straightline), nrow(dat_all)),
  sprintf("- Experience (analytic sample): %s", paste(sprintf("%s %d (%.1f%%)", exp_dist$category, exp_dist$n, exp_dist$pct), collapse = "; ")),
  sprintf("- Frequency (analytic sample): %s", paste(sprintf("%s %d (%.1f%%)", freq_dist$category, freq_dist$n, freq_dist$pct), collapse = "; ")),
  sprintf("- Mardia: skew b1p = %.2f (p %s), kurtosis b2p = %.2f, z = %.2f (p %s); KMO = %.3f", mard$b1p, p_fmt(mard$p.skew), mard$b2p, mard$kurtosis, p_fmt(mard$p.kurt), kmo$MSA),
  sprintf("- CFA (CE, AD): scaled chi2(%d) = %.2f, p = %s; robust CFI = %.3f; robust TLI = %.3f; robust RMSEA = %.3f [%.3f, %.3f]; SRMR = %.3f",
          fm_cfa["df.scaled"], fm_cfa["chisq.scaled"], p_fmt(fm_cfa["pvalue.scaled"]), fm_cfa["cfi.robust"], fm_cfa["tli.robust"],
          fm_cfa["rmsea.robust"], fm_cfa["rmsea.ci.lower.robust"], fm_cfa["rmsea.ci.upper.robust"], fm_cfa["srmr"]),
  sprintf("- Loadings: %s", paste(sprintf("%s %.3f", load$rhs, load$est.std), collapse = "; ")),
  sprintf("- Reliability: CE alpha = %.3f (Spearman-Brown %.3f, inter-item r %.3f), omega = %.3f, AVE = %.3f; AD alpha = %.3f (Spearman-Brown %.3f, inter-item r %.3f), omega = %.3f, AVE = %.3f",
          alpha_ce, sbr(r89), r89, omega["CE"], ave["CE"], alpha_ad, sbr(r67), r67, omega["AD"], ave["AD"]),
  sprintf("- With Q10 as a third dependency indicator: Q10 loading = %.3f; AVE CE = %.3f, AD = %.3f", q10_load, ave3["CE"], ave3["AD"]),
  sprintf("- Latent r(CE, AD) = %.3f [%.3f, %.3f]; HTMT = %.3f; r = 1 test: chi2 diff(1) = %.2f, p = %s",
          r_ce_ad$est.std, r_ce_ad$ci.lower, r_ce_ad$ci.upper, htmt["AD", "CE"], lrt_one[2, "Chisq diff"], p_fmt(lrt_one[2, "Pr(>Chisq)"])),
  sprintf("- CMV: Harman first component = %.1f%%; ULMC method loading = %.3f (%.1f%% of item variance), robust CFI change = %.3f; CE-AD r %.3f -> %.3f [%.3f, %.3f] with ULMC factor",
          harman_pct, m_load[1], 100 * mean(m_load^2), fmt_m["cfi.robust"] - fmt_t["cfi.robust"], r_t$est.std, r_m$est.std, r_m$ci.lower, r_m$ci.upper),
  sprintf("- CMV marker approaches (marker = AI-induced remorse): Lindell-Whitney r_M = %.3f, adjusted r(CE, AD) = %.3f (p %s); CFA marker technique r(CE, AD) = %.3f [%.3f, %.3f]",
          r_mk, r_lw, p_fmt(p_lw), r_cmk$est.std, r_cmk$ci.lower, r_cmk$ci.upper),
  sprintf("- Lindell-Whitney applied to all correlations (scale scores, FREQ, EXP): scale-score beta(CE) %.3f -> %.3f [%.3f, %.3f], p %s",
          lw_b0$est.std, lw_b1$est.std, lw_b1$ci.lower, lw_b1$ci.upper, p_fmt(lw_b1$pvalue)),
  sprintf("- Conservative Lindell-Whitney (r_M = %.3f): adjusted r = %.3f (p %s); scale-score beta(CE) = %.3f [%.3f, %.3f] (p %s)",
          r_mk2, r_lw2, p_fmt(p_lw2), lw_b2$est.std, lw_b2$ci.lower, lw_b2$ci.upper, p_fmt(lw_b2$pvalue)),
  sprintf("- H1 model fit: scaled chi2(%d) = %.2f, p = %s; robust CFI = %.3f; robust TLI = %.3f; robust RMSEA = %.3f [%.3f, %.3f]; SRMR = %.3f",
          fm_h1["df.scaled"], fm_h1["chisq.scaled"], p_fmt(fm_h1["pvalue.scaled"]), fm_h1["cfi.robust"], fm_h1["tli.robust"],
          fm_h1["rmsea.robust"], fm_h1["rmsea.ci.lower.robust"], fm_h1["rmsea.ci.upper.robust"], fm_h1["srmr"]),
  sprintf("- H1: CE -> AD  B = %.3f, robust SE = %.3f, z = %.2f, p = %s, beta = %.3f [%.3f, %.3f]; BC bootstrap 95%% CI for B [%.3f, %.3f]",
          h1_ce$est, h1_ce$se, h1_ce$z, h1_ce$pvalue, h1_ce$est.std, h1_ce$ci.lower, h1_ce$ci.upper, h1_ce$bc_ci_lower, h1_ce$bc_ci_upper),
  sprintf("- Single indicator of CE: omega = %.3f; fixed error variance = %.3f", omega["CE"], ev_ce),
  sprintf("- Bootstrap: %d of 5,000 draws non-admissible (%s), retained as in lavaan; BC CI for CE -> AD without them [%.3f, %.3f]",
          n_nonadm, if (n_nonadm == 0) "none" else paste(names(neg_var)[neg_var > 0], "negative in", neg_var[neg_var > 0], "draws", collapse = "; "),
          h1_ce$bc_ci_lower_admissible, h1_ce$bc_ci_upper_admissible),
  sprintf("- Two-step bootstrap (omega re-estimated in each resample): %d usable resamples; BC CI for B [%.3f, %.3f]", sum(b2_ok), bc2[1], bc2[2]),
  sprintf("- R4a (two-indicator latent CE) bootstrap: %d of 5,000 draws non-admissible (%.1f%%)", n_nonadm_r4, 100 * n_nonadm_r4 / 5000),
  sprintf("- Controls: FREQ beta = %.3f (p %s); EXP beta = %.3f (p %s)", h1_fr$est.std, h1_fr$pvalue, h1_ex$est.std, h1_ex$pvalue),
  sprintf("- R2(AD) = %.3f; R2 without CE = %.3f; f2 for CE = %.3f; chi2 diff for CE path = %.2f, p = %s",
          r2_h1, r2_h1_0, f2_ce, lrt_ce[2, "Chisq diff"], p_fmt(lrt_ce[2, "Pr(>Chisq)"])),
  sprintf("- Minimum detectable r at 80%% power: %.3f", r_min),
  sprintf("- Reliability sensitivity (omega .50 to .80): beta %.3f to %.3f; B %.3f to %.3f; largest p = %.2e",
          min(sens_rel$beta), max(sens_rel$beta), min(sens_rel$B), max(sens_rel$B), max(sens_rel$p)),
  sprintf("- Exploratory (Online Resource 2, Table S8): CE -> PB beta = %.3f; indirect beta = %.3f [%.3f, %.3f], BC bootstrap CI for B [%.3f, %.3f]; direct beta = %.3f; indirect share = %.3f; reverse-order indirect beta = %.3f [%.3f, %.3f]",
          me$est.std[1], me$est.std[4], me$ci.lower[4], me$ci.upper[4], me$bc_lower[4], me$bc_upper[4], me$est.std[3], me$est[6],
          re$est.std[4], re$ci.lower[4], re$ci.upper[4]), "",
  "## Robustness", "",
  paste0("- ", robust$check, " [", robust$predictor, "]: beta = ", robust$beta, " [", robust$beta_ci_lower, ", ", robust$beta_ci_upper,
         "], p = ", robust$p, ", n = ", robust$n, ", R2 = ", robust$R2, ", CFI = ", robust$CFI, ", RMSEA = ", robust$RMSEA, ", SRMR = ", robust$SRMR)
)
writeLines(summ, file.path(out_dir, "results_summary.md"))
cat(summ, sep = "\n")

# -------------------------------------------------------------------------------------
section("11. Values as printed in the manuscript")
# Every value is rounded once from the unrounded estimate (rounding the 3-decimal tables again can be off by .01).
fci <- function(e, lo, hi) sprintf("%s [%s, %s]", d2(e), d2(lo), d2(hi))
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
cor_lines <- unlist(lapply(2:length(cvars), function(i) sapply(1:(i - 1), function(j)
  sprintf("- r(%s, %s) = %s%s", cvars[i], cvars[j], d2(ct$r[i, j]), stars(ct$p[i, j])))))
r_ce_pb <- s2[s2$op == "~~" & ((s2$lhs == "CE" & s2$rhs == "PB") | (s2$lhs == "PB" & s2$rhs == "CE")), "est.std"]
b3s   <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))                  # fit indices: three decimals, no leading zero
p_apa <- function(p) ifelse(p < .001, "< .001", paste("=", sub("^0", "", sprintf("%.3f", p))))
fit_apa <- function(f) sprintf("chi2_YB(%d) = %.2f, p %s, robust CFI = %s, robust TLI = %s, robust RMSEA = %s [%s, %s], SRMR = %s",
                               as.integer(round(f[["df.scaled"]])), f[["chisq.scaled"]], p_apa(f[["pvalue.scaled"]]), b3s(f[["cfi.robust"]]),
                               b3s(f[["tli.robust"]]), b3s(f[["rmsea.robust"]]), b3s(f[["rmsea.ci.lower.robust"]]),
                               b3s(f[["rmsea.ci.upper.robust"]]), b3s(f[["srmr"]]))
mv <- c("# Values as printed in the manuscript (generated by analysis/sem_analysis.R)", "",
  "Each value is rounded once from the unrounded estimate. Values bounded by 1 have no leading zero (APA style).", "",
  "## Participants and analytic sample", "",
  sprintf("- Responses: 557; complete: %d; straight-lining respondents excluded: %d (%.1f%% of the complete responses); analytic N = %d",
          nrow(dat_all), sum(dat_all$straightline), 100 * mean(dat_all$straightline), nrow(dat)),
  sprintf("- Experience (analytic sample): %s", paste(sprintf("%s %.1f%% (n = %d)", exp_dist$category, exp_dist$pct, exp_dist$n), collapse = "; ")),
  sprintf("- Frequency (analytic sample): %s", paste(sprintf("%s %.1f%% (n = %d)", freq_dist$category, freq_dist$pct, freq_dist$n), collapse = "; ")),
  sprintf("- CGPA available for %d of %d", sum(!is.na(dat$CGPA)), nrow(dat)),
  sprintf("- Mardia (four focal items): skewness b1p = %.2f, chi2 = %.2f, p %s; kurtosis b2p = %.2f, z = %.2f, p %s",
          mard$b1p, mard$skew, p_apa(mard$p.skew), mard$b2p, mard$kurtosis, p_apa(mard$p.kurt)),
  sprintf("- Smallest correlation detectable with 80%% power (two-tailed alpha = .05): r = %s", d2(r_min)),
  "", "## Table 2: measurement model", "",
  sprintf("- CFA fit: %s", fit_apa(fm_cfa)),
  sprintf("- %s: M (SD) = %.2f (%.2f); standardised loading = %s", load_u$rhs, d[load_u$rhs, "mean"], d[load_u$rhs, "sd"], d2(load_u$est.std)),
  sprintf("- Critical engagement: alpha = %s; Spearman-Brown = %s; inter-item r = %s; omega = %s; AVE = %s; sqrt(AVE) = %s",
          d2(alpha_ce), d2(sbr(r89)), d2(r89), d2(omega["CE"]), d2(ave["CE"]), d2(sqrt(ave["CE"]))),
  sprintf("- AI dependency: alpha = %s; Spearman-Brown = %s; inter-item r = %s; omega = %s; AVE = %s; sqrt(AVE) = %s",
          d2(alpha_ad), d2(sbr(r67)), d2(r67), d2(omega["AD"]), d2(ave["AD"]), d2(sqrt(ave["AD"]))),
  sprintf("- AVE unrounded: CE = %.4f; AD = %.4f", ave["CE"], ave["AD"]),
  sprintf("- Q10 (not in the scale): M (SD) = %.2f (%.2f); loading as a third dependency indicator = %s; AVE of the three-item scale = %s",
          d[q10, "mean"], d[q10, "sd"], d2(q10_load), d2(ave3["AD"])),
  sprintf("- Latent r(CE, AD) = %s; HTMT = %.2f", fci(r_ce_ad$est.std, r_ce_ad$ci.lower, r_ce_ad$ci.upper), htmt["AD", "CE"]),
  sprintf("- Perceived benefits (five items, covariate in R2): alpha = %s", d2(alpha_pb)),
  sprintf("- r = 1 test: Delta chi2(1) = %.2f, p %s; r < .90 test: Delta chi2(1) = %.2f, p %s",
          lrt_one[2, "Chisq diff"], p_apa(lrt_one[2, "Pr(>Chisq)"]), dv[1, "Chisq diff"], p_apa(dv[1, "Pr(>Chisq)"])),
  "", "## Table 3: descriptive statistics and correlations", "",
  sprintf("- %s: M = %.2f, SD = %.2f, n = %d", cvars, d[cvars, "mean"], d[cvars, "sd"], d[cvars, "n"]),
  cor_lines,
  sprintf("- p(EXP, AD_mean) = %s (unadjusted)", sub("^0", "", sprintf("%.3f", ct$p["EXP", "AD_mean"]))),   # lower triangle = unadjusted
  sprintf("- Item skewness %.2f to %.2f; kurtosis %.2f to %.2f (four focal items)", min(d[c(ce_items, ad_items), "skew"]), max(d[c(ce_items, ad_items), "skew"]),
          min(d[c(ce_items, ad_items), "kurtosis"]), max(d[c(ce_items, ad_items), "kurtosis"])),
  "", "## Common-method variance", "",
  sprintf("- Harman: first unrotated component = %.1f%% of the variance of %d items", harman_pct, length(harman_items)),
  sprintf("- Single-factor CFA of the four focal items: robust CFI = %s, robust RMSEA = %s", b3s(fitMeasures(fit_1f, "cfi.robust")), b3s(fitMeasures(fit_1f, "rmsea.robust"))),
  sprintf("- ULMC: CE-AD r = %s without and %s with the method factor; method factor share of item variance = %.0f%%; robust CFI change = %s",
          d2(r_t$est.std), fci(r_m$est.std, r_m$ci.lower, r_m$ci.upper), 100 * mean(m_load^2), b3s(fmt_m["cfi.robust"] - fmt_t["cfi.robust"])),
  sprintf("- Lindell-Whitney: marker r = %s (the smaller of %s with AD and %s with CE); r(CE, AD) %s -> %s, t(%d) = %.2f, p %s",
          d2(r_mk), d2(r_mk_ad), d2(r_mk_ce), d2(r_obs), d2(r_lw), nrow(dat) - 3, t_lw, p_apa(p_lw)),
  sprintf("- Lindell-Whitney applied to all correlations among the scale scores and controls: scale-score beta(CE) %s -> %s, p %s (%.0f%% smaller); R2 %s -> %s",
          d2(lw_b0$est.std), fci(lw_b1$est.std, lw_b1$ci.lower, lw_b1$ci.upper), p_apa(lw_b1$pvalue),
          100 * (1 - lw_b1$est.std / lw_b0$est.std), d2(lw_r2[1]), d2(lw_r2[2])),
  sprintf("- Conservative Lindell-Whitney (r_M = %s, the larger marker correlation): r(CE, AD) %s -> %s, t(%d) = %.2f, p %s; scale-score beta(CE) %s -> %s, p %s (%.0f%% smaller)",
          d2(r_mk2), d2(r_obs), d2(r_lw2), nrow(dat) - 3, t_lw2, p_apa(p_lw2), d2(lw_b0$est.std),
          fci(lw_b2$est.std, lw_b2$ci.lower, lw_b2$ci.upper), p_apa(lw_b2$pvalue), 100 * (1 - lw_b2$est.std / lw_b0$est.std)),
  sprintf("- ULMC: reduction of the latent correlation = %.0f%%", 100 * (1 - r_m$est.std / r_t$est.std)),
  sprintf("- CFA marker: method share of focal-item variance = %.1f%%; r(CE, AD) = %s", 100 * mean(mm_load^2),
          fci(r_cmk$est.std, r_cmk$ci.lower, r_cmk$ci.upper)),
  "", "## Table 4: H1 structural model", "",
  sprintf("- Fit: %s", fit_apa(fm_h1)),
  sprintf("- Single indicator of critical engagement: omega = %s (%.3f); variance of the mean = %.3f; fixed error variance = (1 - %.3f) x %.3f = %.3f",
          d2(omega["CE"]), omega["CE"], var(dat$CE_mean), omega["CE"], var(dat$CE_mean), ev_ce),
  sprintf("- %s: B = %s, SE = %s, z = %s, p = %s, beta = %s, BC bootstrap CI for B = %s (admissible draws only: %s)",
          h1_apa$rhs, h1_apa$B, h1_apa$SE, h1_apa$z, h1_apa$p, h1_apa$beta, h1_apa$bc_ci, h1_apa$bc_ci_admissible),
  sprintf("- R2 = %s; R2 with the CE path fixed to 0 = %s; f2 = %.2f; Delta chi2(1) for the CE path = %.2f, p %s; non-admissible bootstrap draws = %d (%.1f%%)",
          d2(r2_h1), d2(r2_h1_0), f2_ce, lrt_ce[2, "Chisq diff"], p_apa(lrt_ce[2, "Pr(>Chisq)"]), n_nonadm, 100 * n_nonadm / 5000),
  sprintf("- Two-step bootstrap (omega re-estimated in each resample): BC CI for B = [%.2f, %.2f]; usable resamples = %d of 5000",
          bc2[1], bc2[2], sum(b2_ok)),
  sprintf("- Predictor correlation %s-%s = %s", covs$lhs, covs$rhs, d2(covs$est.std)),
  sprintf("- R4a (two-indicator latent critical engagement): non-admissible bootstrap draws = %d of 5000 (%.1f%%)", n_nonadm_r4, 100 * n_nonadm_r4 / 5000),
  "", "## Reliability sensitivity (Online Resource 2, Table S9 and Fig. S1)", "",
  sprintf("- Assumed omega = %.2f: error variance = %.3f; B = %.2f, SE = %.2f, z = %.2f, p %s; beta = %s; R2 = %s",
          sens_rel$reliability, sens_rel$error_variance, sens_rel$B, sens_rel$SE, sens_rel$z, p_apa(sens_rel$p),
          fci(sens_rel$beta, sens_rel$beta_ci_lower, sens_rel$beta_ci_upper), d2(sens_rel$R2)),
  sprintf("- Range across the grid: beta %s to %s; B %.2f to %.2f; SE %.2f to %.2f; all p %s; omega in the two-step bootstrap (2.5%% to 97.5%%): %s to %s",
          d2(min(sens_rel$beta)), d2(max(sens_rel$beta)), min(sens_rel$B), max(sens_rel$B), min(sens_rel$SE), max(sens_rel$SE),
          p_apa(max(sens_rel$p)), d2(om_rng[1]), d2(om_rng[2])),
  "", "## Table 5: robustness checks", "",
  sprintf("- %s [%s]: n = %d, beta = %s, p = %s, R2 = %s", robust$check, robust$predictor, robust$n, robust$beta_apa, robust$p, robust$R2_apa),
  sprintf("- R2 model: latent r(CE, PB) = %s", d2(r_ce_pb)),
  sprintf("- R4c: Spearman-Brown coefficient used = %s; R6: omega in the full sample = %s", d2(sbr(r89)), d2(omega_all[["CE"]])),
  "", "## Exploratory analysis (Online Resource 2, Table S8)", "",
  sprintf("- %s: B = %.2f, SE = %.2f, p = %s, beta = %s, BC bootstrap CI for B = [%.2f, %.2f]",
          c("a (CE -> PB)", "b (PB -> AD)", "c' (CE -> AD, direct)", "indirect (a x b)", "total (c' + a x b)"),
          me$est[1:5], me$se[1:5], p_fmt(me$pvalue[1:5]), fci(me$est.std[1:5], me$ci.lower[1:5], me$ci.upper[1:5]), me$bc_lower[1:5], me$bc_upper[1:5]),
  sprintf("- Indirect share of the total association = %.0f%%; R2: PB = %s, AD = %s; non-admissible bootstrap draws = %d",
          100 * me$est[6], d2(r2_med["PB"]), d2(r2_med["AD"]), length(nonadm_med)),
  sprintf("- Reverse ordering (PB -> CE -> AD), identical fit: %s: beta = %s, p = %s",
          c("a2 (PB -> CE)", "b2 (CE -> AD)", "c2' (PB -> AD, direct)", "indirect", "total"),
          fci(re$est.std, re$ci.lower, re$ci.upper), p_fmt(re$pvalue)))
writeLines(mv, file.path(out_dir, "manuscript_values.md"))
cat(mv, sep = "\n")

# -------------------------------------------------------------------------------------
section("12. Online Resource 2: item-level statistics, full model output and the exploratory analysis")
# Writes online_resource_2_tables.tex, typeset by manuscript/online_resource_2.tex (uploaded as ESM_2).
# As in Section 11, every value is rounded once from the unrounded estimate.
tm   <- function(s) gsub("-", "$-$", s, fixed = TRUE)                    # typographic minus (numbers only)
n2   <- function(x) ifelse(is.na(x), "--", tm(sub("^-(0\\.0+)$", "\\1", sprintf("%.2f", x))))   # statistics that can exceed 1 (no "-0.00")
b2   <- function(x) ifelse(is.na(x), "--", tm(d2(x)))                    # statistics bounded by 1: no leading zero
b3   <- function(x) ifelse(is.na(x), "--", tm(sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))))   # fit indices
p3   <- function(p) ifelse(is.na(p), "--", ifelse(p < .001, "$<.001$", sub("^0", "", sprintf("%.3f", p))))
p_eq <- function(p) ifelse(p < .001, "$p < .001$", sprintf("$p = %s$", sub("^0", "", sprintf("%.3f", p))))
ci_b <- function(e, lo, hi, f = b2) sprintf("%s [%s, %s]", f(e), f(lo), f(hi))
tr   <- function(...) paste0(paste(..., sep = " & "), " \\\\")
grp  <- function(text, ncol) sprintf("\\multicolumn{%d}{l}{\\emph{%s}} \\\\", ncol, text)
table_tex <- function(env, caption, label, colspec, header, body, note, size = "\\footnotesize") {
  c(sprintf("\\begin{%s}%s", env, if (env == "table") "[H]" else ""), "\\centering", sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label), size, "\\begin{threeparttable}", sprintf("\\begin{tabular}{%s}", colspec), "\\toprule",
    header, "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}\\footnotesize", paste("\\item", note),
    "\\end{tablenotes}", "\\end{threeparttable}", sprintf("\\end{%s}", env), "")
}
fit_idx <- function(fit) {                       # robust indices for MLR, scaled indices for WLSMV
  nm <- fm_names[-4]
  if (lavInspect(fit, "options")$estimator == "DWLS") nm <- sub("robust", "scaled", nm)
  unname(fitMeasures(fit, nm))
}
chi_df <- function(f) sprintf("%s (%d)", n2(f[1]), as.integer(round(f[2])))

esm_items <- data.frame(
  var = c(ce_items, ad_items, q10, "FREQ", "EXP", pb_items, ar_items[c(4, 5, 1, 2, 3)], "CGPA"),
  label = c("Q8 Starting point for own ideas", "Q9 Cross-checking with other sources", "Q6 Reliance for decisions",
            "Q7 AI instead of own problem-solving", "Q10 Acceptance without verification", "Q5 Frequency of use (1--5)",
            "Q4 AI-use experience (1--4)", "Q17 Saves time", "Q18 Confidence in problem-solving", "Q19 Learning efficiency",
            "Q20 Satisfaction with learning", "Q21 Collaborative learning tool", "Q11 Fear of weakened thinking",
            "Q12 Reduced recall", "Q14 Remorse for avoiding challenges", "Q15 Regret over dependency",
            "Q16 Regret after tasks", "CGPA (0--4)"),
  group = c(rep("Critical engagement", 2), rep("AI dependency", 2), "Not in the dependency scale (R3, R9 and R10)", rep("Control variables", 2),
            rep("Perceived benefits (covariate in R2)", 5), rep("AI-induced remorse (marker variable)", 5),
            "Academic performance (covariate in R8)"),
  k = c(rep(7, 5), 5, 4, rep(7, 10), NA), stringsAsFactors = FALSE)
item_lab <- setNames(esm_items$label, esm_items$var)

# Table S1: item descriptives and response distributions
de <- psych::describe(dat[, esm_items$var])
s1 <- character(); g_prev <- ""
for (i in seq_len(nrow(esm_items))) {
  v <- esm_items$var[i]; x <- dat[[v]]; x <- x[!is.na(x)]; k <- esm_items$k[i]
  if (esm_items$group[i] != g_prev) { s1 <- c(s1, grp(esm_items$group[i], 13)); g_prev <- esm_items$group[i] }
  if (!is.na(k)) stopifnot(all(x %in% seq_len(k)))
  pct <- sapply(1:7, function(j) if (is.na(k) || j > k) "--" else sprintf("%.1f", 100 * mean(x == j)))
  s1 <- c(s1, tr(esm_items$label[i], length(x), sprintf("%.2f", mean(x)), sprintf("%.2f", sd(x)),
                 n2(de[v, "skew"]), n2(de[v, "kurtosis"]), paste(pct, collapse = " & ")))
}
cat("Item descriptives (Table S1):\n")
print(data.frame(item = esm_items$label, n = de$n, M = r3(de$mean), SD = r3(de$sd), skew = r3(de$skew), kurtosis = r3(de$kurtosis)), row.names = FALSE)

# Table S2: correlation matrix of all items and covariates (Pearson; lower triangle of corr.test p is unadjusted)
ct_i <- psych::corr.test(dat[, esm_items$var], use = "pairwise")
nv <- nrow(esm_items)
s2 <- sapply(seq_len(nv), function(i) tr(sprintf("%d. %s", i, esm_items$label[i]), paste(sapply(seq_len(nv - 1), function(j)
  if (j < i) paste0(b2(ct_i$r[i, j]), stars(ct_i$p[i, j])) else ""), collapse = " & ")))

# Table S3: measurement model in full
fac_label <- c(CE = "Critical engagement (two-factor CFA, Table 2 of the article)", AD = "AI dependency (two-factor CFA)",
               PB = "Perceived benefits (measurement part of robustness model R2)")
meas_rows <- function(fit, factors) {
  pe <- parameterEstimates(fit); ss <- standardizedSolution(fit); rsq <- lavInspect(fit, "rsquare")
  unlist(lapply(factors, function(f) c(grp(fac_label[[f]], 8), sapply(pe$rhs[pe$op == "=~" & pe$lhs == f], function(it) {
    e <- pe[pe$op == "=~" & pe$lhs == f & pe$rhs == it, ]; s <- ss[ss$op == "=~" & ss$lhs == f & ss$rhs == it, ]
    fixed <- is.na(e$z)                          # first loading fixed to 1 (marker scaling)
    tr(item_lab[[it]], if (fixed) "1.00\\tnote{a}" else n2(e$est), if (fixed) "--" else n2(e$se), if (fixed) "--" else n2(e$z),
       if (fixed) "--" else p3(e$pvalue), ci_b(s$est.std, s$ci.lower, s$ci.upper),
       n2(pe$est[pe$op == "~~" & pe$lhs == it & pe$rhs == it]), b2(rsq[[it]]))
  }))))
}
s3 <- c(meas_rows(fit_cfa, c("CE", "AD")), meas_rows(fit_r2, "PB"))

# Tables S4 and S5: every regression path and the fit of the primary and robustness models
mod_list <- list(
  list("Primary model (H1)", fit_h1), list("R1 No control variables", fit_r1),
  list("R2 Adds perceived benefits (latent covariate)", fit_r2),
  list("R3 Four-item dependency scale including the frequency item", fit_r3),
  list("R4 Critical engagement as a latent factor (two indicators)", fit_r4a),
  list("R4 Both constructs as reliability-corrected single indicators", fit_r4b),
  list("R4 Correction based on the Spearman--Brown coefficient", fit_r4c),
  list("R5 Items Q8 and Q9 entered together", fit_r5), list("R5 Item Q9 entered alone", fit_r5c),
  list("R6 All respondents, including straight-lining respondents", fit_r6), list("R7 Ordinal estimation (WLSMV)", fit_r7),
  list("R8 Adds self-reported CGPA", fit_r8), list("R9 Dependency scale including Q10", fit_r9),
  list("R10 Outcome: acceptance without verification (Q10)", fit_r10))
pred_lab <- c(CE = "Critical engagement", CEs = "Critical engagement (single indicator)", FREQ = "Frequency of use",
              EXP = "AI-use experience", PB = "Perceived benefits", CGPA = "CGPA",
              CE_Q8_starting_point = "Q8 Starting point for own ideas", CE_Q9_crosscheck = "Q9 Cross-checking with other sources")
s4 <- unlist(lapply(mod_list, function(m) {
  pe <- parameterEstimates(m[[2]]); ss <- standardizedSolution(m[[2]]); pr <- pe[pe$op == "~", ]
  c(grp(m[[1]], 6), sapply(seq_len(nrow(pr)), function(j) {
    s <- ss[ss$op == "~" & ss$lhs == pr$lhs[j] & ss$rhs == pr$rhs[j], ]
    tr(pred_lab[[pr$rhs[j]]], n2(pr$est[j]), n2(pr$se[j]), n2(pr$z[j]), p3(pr$pvalue[j]), ci_b(s$est.std, s$ci.lower, s$ci.upper))
  }))
}))
h4 <- "Predictor & $B$ & $SE$ & $z$ & $p$ & $\\beta$ [95\\% CI] \\\\"
s4_tex <- c("\\begin{footnotesize}", "\\begin{longtable}{lccccc}",
  "\\caption{Regression estimates for the primary model and every robustness model}\\label{tab:S4} \\\\",
  "\\toprule", h4, "\\midrule", "\\endfirsthead",
  "\\multicolumn{6}{l}{\\emph{Table S4 (continued)}} \\\\", "\\toprule", h4, "\\midrule", "\\endhead",
  "\\midrule", "\\multicolumn{6}{r}{\\emph{continued on next page}} \\\\", "\\endfoot", "\\bottomrule", "\\endlastfoot",
  s4, "\\end{longtable}", "\\end{footnotesize}",
  paste("{\\footnotesize \\emph{Note.} Outcome: AI dependency, a latent factor measured by Q6 and Q7, except in R3 (four-item factor:",
        "Q5, Q6, Q7 and Q10), the second R4 model (reliability-corrected single indicator), R9 (three-item factor: Q6, Q7 and Q10) and R10",
        "(Q10, observed). Critical engagement is the mean of Q8 and Q9 with its error variance fixed at $(1 - \\omega)$ times its variance,",
        "except in the first R4 model and R7 (latent factor with two indicators) and in R5 (items entered as observed predictors); the",
        "third R4 model fixes the error variance with the Spearman--Brown coefficient. MLR estimates with robust standard errors, except",
        "R7 (WLSMV, items treated as ordinal). Standardised 95\\% CIs by the delta method.",
        sprintf("$n = %d$, except R6 ($n = %d$, all complete responses) and R8 ($n = %d$).", nrow(dat), nrow(dat_all), nrow(d_cg)),
        "The specifications are described in Section 3.4 and Table 5 of the article; their fit indices are in Table S5.\\par}"), "")
fit_line <- function(label, fit) {
  f <- fit_idx(fit); pe <- parameterEstimates(fit); dv <- unique(pe$lhs[pe$op == "~"])
  tr(label, lavInspect(fit, "nobs"), chi_df(f), p3(f[3]), b3(f[4]), b3(f[5]), ci_b(f[6], f[7], f[8], b3), b3(f[9]),
     if (length(dv) == 1) b2(lavInspect(fit, "rsquare")[[dv]]) else "--")
}
s5 <- c(fit_line("Two-factor CFA (Table 2)", fit_cfa), sapply(mod_list, function(m) fit_line(m[[1]], m[[2]])))
cat("Test statistic used by lavaan for WLSMV (R7):", paste(lavInspect(fit_r7, "options")$test, collapse = " + "), "\n")

# Tables S6 and S7: common-method variance diagnostics in full
cmv_line <- function(label, fit, r = NULL) {
  f <- fit_idx(fit)
  tr(label, chi_df(f), p3(f[3]), b3(f[4]), b3(f[5]), ci_b(f[6], f[7], f[8], b3), b3(f[9]),
     if (is.null(r)) "--" else ci_b(r$est.std, r$ci.lower, r$ci.upper))
}
s6 <- c(cmv_line("Two-factor CFA, four focal items (Table 2)", fit_cfa, r_ce_ad),
        cmv_line("Single-factor CFA, four focal items", fit_1f),
        cmv_line("Four-factor CFA, 14 seven-point items", fit_t, r_t),
        cmv_line("Four-factor CFA with ULMC", fit_m, r_m),
        cmv_line("Two-factor CFA with CFA marker factor", fit_mk, r_cmk))
s7 <- c(
  tr(sprintf("Harman: variance explained by the first unrotated component (%d items)", length(harman_items)), sprintf("%.1f\\%%", harman_pct)),
  tr("Single- versus two-factor model of the four focal items: scaled $\\Delta\\chi^2(1)$",
     sprintf("%s, %s", n2(lrt_one[2, "Chisq diff"]), p_eq(lrt_one[2, "Pr(>Chisq)"]))),
  tr("ULMC: standardised method loadings (range)", sprintf("%s to %s", b2(min(m_load)), b2(max(m_load)))),
  tr("ULMC: mean share of item variance due to the method factor", sprintf("%.0f\\%%", 100 * mean(m_load^2))),
  tr("ULMC: change in robust CFI", b3(fmt_m["cfi.robust"] - fmt_t["cfi.robust"])),
  tr("Lindell--Whitney: observed $r$(CE, AD), scale means", b2(r_obs)),
  tr("Lindell--Whitney: marker correlation $r_M$", b2(r_mk)),
  tr("Lindell--Whitney: adjusted $r$(CE, AD)", sprintf("%s, $t(%d) = %s$, %s", b2(r_lw), nrow(dat) - 3, n2(t_lw), p_eq(p_lw))),
  tr("Scale-score regression with both controls: $\\beta$ of critical engagement, unadjusted", ci_b(lw_b0$est.std, lw_b0$ci.lower, lw_b0$ci.upper)),
  tr("The same regression with every correlation adjusted by $r_M$", sprintf("%s, %s", ci_b(lw_b1$est.std, lw_b1$ci.lower, lw_b1$ci.upper), p_eq(lw_b1$pvalue))),
  tr("Conservative variant: $r_M$ = the larger marker correlation", b2(r_mk2)),
  tr("Conservative variant: adjusted $r$(CE, AD)", sprintf("%s, $t(%d) = %s$, %s", b2(r_lw2), nrow(dat) - 3, n2(t_lw2), p_eq(p_lw2))),
  tr("Conservative variant: adjusted regression, $\\beta$ of critical engagement", sprintf("%s, %s", ci_b(lw_b2$est.std, lw_b2$ci.lower, lw_b2$ci.upper), p_eq(lw_b2$pvalue))),
  tr("CFA marker: standardised method loadings on the focal items (range)", sprintf("%s to %s", b2(min(mm_load)), b2(max(mm_load)))),
  tr("CFA marker: mean share of focal-item variance due to the marker factor", sprintf("%.0f\\%%", 100 * mean(mm_load^2))))

# Table S8: exploratory analysis (Section 9)
lab8 <- c(a = "Critical engagement $\\rightarrow$ perceived benefits ($a$)", b = "Perceived benefits $\\rightarrow$ AI dependency ($b$)",
          cp = "Critical engagement $\\rightarrow$ AI dependency, direct ($c'$)", ind = "Indirect association ($a \\times b$)",
          tot = "Total association ($c' + a \\times b$)", a2 = "Perceived benefits $\\rightarrow$ critical engagement",
          b2 = "Critical engagement $\\rightarrow$ AI dependency", cp2 = "Perceived benefits $\\rightarrow$ AI dependency, direct",
          ind2 = "Indirect association", tot2 = "Total association")
row8 <- function(r, boot) tr(lab8[[r$param]], n2(r$est), n2(r$se), n2(r$z), p3(r$pvalue), ci_b(r$est.std, r$ci.lower, r$ci.upper),
                             if (boot) sprintf("[%s, %s]", n2(r$bc_lower), n2(r$bc_upper)) else "--")
s8 <- c(grp("A. Perceived benefits as an intervening variable", 7), sapply(c("a", "b", "cp", "ind", "tot"), function(p) row8(me[me$param == p, ], TRUE)),
        grp("B. Equivalent model with the reverse ordering", 7), sapply(c("a2", "b2", "cp2", "ind2", "tot2"), function(p) row8(re[re$param == p, ], FALSE)))
boot_note <- if (length(nonadm_med) == 0) "none of the draws was non-admissible" else
  sprintf("%d non-admissible draws were retained; without them, the CI for the indirect association is [%s, %s]",
          length(nonadm_med), n2(ind_bc_adm[1]), n2(ind_bc_adm[2]))

# Table S9 and Fig. S1: sensitivity of the primary estimate to the reliability assumed for critical engagement (Section 8)
s9 <- sapply(seq_len(nrow(sens_rel)), function(i) with(sens_rel[i, ], tr(sprintf("%.2f", reliability), sprintf("%.3f", error_variance),
  n2(B), n2(SE), n2(z), p3(p), ci_b(beta, beta_ci_lower, beta_ci_upper), b2(R2))))
s9_tex <- c(
  table_tex("table", "Sensitivity of the primary estimate to the reliability assumed for critical engagement", "tab:S9", "cccccccc",
    "Assumed $\\omega$ & Error variance & $B$ & $SE$ & $z$ & $p$ & $\\beta$ [95\\% CI] & $R^2$ \\\\", s9,
    sprintf(paste("\\emph{Note.} $N = %d$. Each row re-estimates the primary model with the error variance of critical engagement",
                  "fixed at $(1 - \\omega) \\times %s$, where %s is the variance of the mean of Q8 and Q9. The estimated reliability is",
                  "$\\omega = %s$. In the two-step bootstrap, 95\\%% of the $\\omega$ estimates lay between %s and %s (Online Resource 1,",
                  "Section 6). MLR estimates with robust standard errors; standardised 95\\%% CIs by the delta method. $R^2$ refers to",
                  "AI dependency. Fig.~\\ref{fig:S1} plots the estimates."),
            nrow(dat), sprintf("%.3f", var(dat$CE_mean)), sprintf("%.3f", var(dat$CE_mean)), sub("^0", "", sprintf("%.3f", omega[["CE"]])),
            b2(om_rng[1]), b2(om_rng[2]))),
  "\\begin{figure}[H]", "\\centering", "\\includegraphics[width=\\textwidth]{../analysis/output/fig_S1_reliability_sensitivity.pdf}",
  sprintf(paste("\\caption{Coefficient of critical engagement in the primary model for assumed reliabilities of .50 to .80, with 95\\%%",
                "confidence intervals: (a) unstandardised, with robust standard errors; (b) standardised. The dotted line marks the",
                "estimated $\\omega$ (%s). The shaded band covers the central 95\\%% of the $\\omega$ estimates in the two-step bootstrap.",
                "The dashed line marks zero}"), b2(omega[["CE"]])),
  "\\label{fig:S1}", "\\end{figure}", "")

esm_tex <- c("% Tables of Online Resource 2, generated by analysis/sem_analysis.R (Section 12). Do not edit by hand.", "",
  table_tex("table", "Item-level descriptive statistics and response distributions", "tab:S1", "lrrrrrrrrrrrr",
    c("& & & & & & \\multicolumn{7}{c}{Responses at each scale point (\\%)} \\\\", "\\cmidrule(l){7-13}",
      "Item & $n$ & $M$ & $SD$ & Skew. & Kurt. & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\\\"), s1,
    paste("\\emph{Note.} Item wording is given in Appendix Table A1 of the article. Seven-point items range from 1 (lowest) to 7",
          "(highest agreement or reliance); Q5: 1 = \\emph{never} to 5 = \\emph{always}; Q4: 1 = less than 6 months to 4 = more than",
          "2 years; CGPA: self-reported, 0--4. Skewness and excess kurtosis as computed by the \\texttt{describe} function of psych.",
          "Percentages may not sum to 100 because of rounding."), size = "\\footnotesize\\setlength{\\tabcolsep}{3pt}"),
  table_tex("sidewaystable", "Correlations among all items and covariates used in the analyses", "tab:S2",
    paste0("l", strrep("c", nv - 1)), tr("Variable", paste(seq_len(nv - 1), collapse = " & ")), s2,
    paste(sprintf("\\emph{Note.} Pearson correlations, $N = %d$ (pairs involving CGPA: $n = %d$). Means, standard deviations and response", nrow(dat), nrow(d_cg)),
          "distributions are given in Table S1. * $p < .05$; ** $p < .01$; *** $p < .001$ (two-tailed, unadjusted)."),
    size = "\\footnotesize\\setlength{\\tabcolsep}{3pt}"),
  "\\clearpage",
  table_tex("table", "Measurement model: full parameter estimates", "tab:S3", "lccccccc",
    "Item & $B$ & $SE$ & $z$ & $p$ & $\\lambda$ [95\\% CI] & $\\theta$ & $R^2$ \\\\", s3,
    c(paste(sprintf("\\emph{Note.} MLR estimates, $N = %d$. $B$ = unstandardised loading with robust standard error; $\\lambda$ = standardised", nrow(dat)),
            "loading with delta-method 95\\% CI; $\\theta$ = unstandardised residual variance; $R^2$ = proportion of item variance",
            "explained by the factor. Fit indices are given in Table S5. In the primary structural model, critical engagement is the mean",
            sprintf("of Q8 and Q9 with its error variance fixed at $(1 - \\omega)$ times its variance: $(1 - %s) \\times %s = %s$.",
                    sub("^0", "", sprintf("%.3f", omega[["CE"]])), sprintf("%.3f", var(dat$CE_mean)), sprintf("%.3f", ev_ce))),
      "[a] Fixed to 1 to set the scale of the factor.")),
  s4_tex,
  table_tex("table", "Fit indices and explained variance for the primary model and every robustness model", "tab:S5",
    ">{\\raggedright\\arraybackslash}p{4.6cm}cccccccc",
    "Model & $n$ & $\\chi^2$ ($df$) & $p$ & CFI & TLI & RMSEA [90\\% CI] & SRMR & $R^2$ \\\\", s5,
    paste("\\emph{Note.} For MLR models, $\\chi^2$ is the Yuan--Bentler scaled statistic and CFI, TLI and RMSEA are the robust versions;",
          "for R7 (WLSMV), the scaled-and-shifted statistic and the scaled indices are reported. $R^2$ refers to the outcome of each",
          "model. Models with $df = 0$ are saturated, so their fit cannot be tested. The RMSEA is unreliable in models with few degrees",
          "of freedom.")),
  table_tex("table", "Common-method variance: fit of the compared models", "tab:S6", ">{\\raggedright\\arraybackslash}p{4.2cm}ccccccc",
    "Model & $\\chi^2$ ($df$) & $p$ & CFI & TLI & RMSEA [90\\% CI] & SRMR & Latent $r$ [95\\% CI] \\\\", s6,
    paste(sprintf("\\emph{Note.} MLR estimates, $N = %d$; Yuan--Bentler scaled $\\chi^2$ and robust fit indices. Latent $r$ is the correlation", nrow(dat)),
          "between critical engagement and AI dependency. The four-factor CFA comprises critical engagement, AI dependency, perceived",
          "benefits and AI-induced remorse (residual covariance of Q11 and Q12 freed). The ULMC loads all 14 items on a method factor",
          "with equal unstandardised loadings, uncorrelated with the trait factors. The CFA marker model adds the five remorse items as",
          "a marker factor, with equal unstandardised loadings on the four focal items and no correlation with the focal factors.")),
  table_tex("table", "Common-method variance: further diagnostics", "tab:S7", "ll", "Diagnostic & Result \\\\", s7,
    paste(sprintf("\\emph{Note.} $N = %d$. Harman's test uses the 14 seven-point items of the measurement models and the frequency item.", nrow(dat)),
          sprintf("The marker variable is AI-induced remorse; $r_M$ is the smaller of its correlations with critical engagement (%s) and AI dependency (%s), computed from scale means.",
                  b2(r_mk_ce), b2(r_mk_ad)),
          "The adjusted regression applies $r_M$ to every correlation among the two scale scores, frequency of use and AI-use experience",
          "and re-estimates the regression of dependency on the other three (Malhotra et al., 2006); scale scores are not corrected for",
          "unreliability, and the CIs are normal-theory. The method loadings are equal in unstandardised form; their standardised values",
          "differ across items.")),
  table_tex("table", "Exploratory analysis: indirect association between critical engagement and AI dependency through perceived benefits",
    "tab:S8", ">{\\raggedright\\arraybackslash}p{6.2cm}cccccc",
    "Path & $B$ & $SE$ & $z$ & $p$ & $\\beta$ [95\\% CI] & Bootstrap 95\\% CI for $B$ \\\\", s8,
    paste0("\\emph{Note.} Exploratory analysis, not a hypothesis test; $N = ", nrow(dat), "$. Critical engagement is the reliability-corrected ",
           "single indicator of the primary model. Frequency of use and AI-use experience are controlled in ",
           "every equation. Both models re-parameterise robustness model R2 and therefore have its fit (Table S5). MLR estimates with robust ",
           "standard errors; standardised 95\\% CIs by the delta method. Bootstrap CIs are bias-corrected (5,000 resamples, maximum likelihood; ",
           boot_note, ") and were computed for Panel A only. ",
           sprintf("In Panel A, the indirect association accounts for %.0f\\%% of the total association; $R^2$ = %s for perceived benefits and %s for AI dependency. ",
                   100 * me$est[me$param == "prop"], b2(r2_med[["PB"]]), b2(r2_med[["AD"]])),
           "Because the two orderings fit the data equally well, the indirect association shows how much of the association is shared with ",
           "perceived benefits; it does not establish mediation.")),
  s9_tex)
writeLines(esm_tex, file.path(out_dir, "online_resource_2_tables.tex"))
cat("Online Resource 2 tables written to", file.path(out_dir, "online_resource_2_tables.tex"), "\n")
sink()
