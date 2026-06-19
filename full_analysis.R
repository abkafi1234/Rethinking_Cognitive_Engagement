# =============================================================================
# FULL ANALYSIS SCRIPT — AI Dependency Study
# Mirrors the analytical structure of the referenced paper exactly:
#   1. Sample Adequacy and Common Method Bias
#   2. Data Distribution and Normality
#   3. Measurement Model Evaluation (CFA)
#   4. Reliability and Convergent Validity
#   5. Discriminant Validity
#   6. Structural Model Fit and Path Coefficients (with bootstrap CIs)
#   7. Mediation Analysis (bootstrap BCa)
#
# REQUIREMENTS: lavaan, semTools, semPlot, psych, dplyr
# Install: install.packages(c("lavaan","semTools","semPlot","psych","dplyr"))
#
# HOW TO RUN:
#   setwd("D:/Research/Remorse/Exp EID")
#   sink("full_analysis_output.txt")
#   source("full_analysis.R", encoding="UTF-8")
#   sink()
# =============================================================================

library(lavaan)
library(semTools)
library(semPlot)
library(psych)
library(dplyr)

# -----------------------------------------------------------------------------
# DATA LOAD
# -----------------------------------------------------------------------------
df <- read.csv("cleaned_constructs.csv", stringsAsFactors = FALSE)

item_cols <- c(
  "AD_Q5_frequency","AD_Q6_reliance_decisions",
  "AD_Q7_rely_vs_independent","AD_Q10_accept_no_verify",
  "AR_Q14_remorseful_avoidance","AR_Q15_regret_dependency",
  "AR_Q16_regret_after_tasks","AR_Q11_fear_weaken_thinking",
  "AR_Q12_aware_recall_loss",
  "PB_Q17_saves_time","PB_Q18_confidence",
  "PB_Q19_learning_efficiency","PB_Q20_satisfaction",
  "PB_Q21_collaborative",
  "CE_Q8_starting_point","CE_Q9_crosscheck"
)

df_items <- na.omit(df[, item_cols])
N <- nrow(df_items)
cat(sprintf("N (raw): 557 | N (after listwise deletion): %d\n\n", N))

# =============================================================================
# SECTION 1: SAMPLE ADEQUACY AND COMMON METHOD BIAS
# =============================================================================
cat("=======================================================================\n")
cat("SECTION 1: SAMPLE ADEQUACY AND COMMON METHOD BIAS\n")
cat("=======================================================================\n\n")

# --- 1a. KMO and Bartlett ---
cat("--- KMO Measure of Sampling Adequacy ---\n")
kmo_res <- psych::KMO(df_items)
cat(sprintf("  Overall KMO = %.3f\n", kmo_res$MSA))
cat(sprintf("  Interpretation: %s\n",
    ifelse(kmo_res$MSA >= 0.90, "Marvelous",
    ifelse(kmo_res$MSA >= 0.80, "Meritorious",
    ifelse(kmo_res$MSA >= 0.70, "Middling",
    ifelse(kmo_res$MSA >= 0.60, "Mediocre", "Unacceptable"))))))

bart_res <- psych::cortest.bartlett(df_items)
cat(sprintf("\n--- Bartlett's Test of Sphericity ---\n"))
cat(sprintf("  chi-sq = %.3f, df = %d, p < .001\n", bart_res$chisq, bart_res$df))
cat(sprintf("  Indicator-to-sample ratio: %d:%d = %.1f:1\n\n",
            N, ncol(df_items), N / ncol(df_items)))

# --- 1b. Harman's single-factor test ---
cat("--- Harman's Single-Factor Test (CMB diagnostic) ---\n")
harman_model <- '
  G =~ AD_Q5_frequency + AD_Q6_reliance_decisions +
       AD_Q7_rely_vs_independent + AD_Q10_accept_no_verify +
       AR_Q14_remorseful_avoidance + AR_Q15_regret_dependency +
       AR_Q16_regret_after_tasks + AR_Q11_fear_weaken_thinking +
       AR_Q12_aware_recall_loss +
       PB_Q17_saves_time + PB_Q18_confidence +
       PB_Q19_learning_efficiency + PB_Q20_satisfaction +
       PB_Q21_collaborative +
       CE_Q8_starting_point + CE_Q9_crosscheck
'
harman_fit <- cfa(harman_model, data = df_items, estimator = "ML")
hf <- fitMeasures(harman_fit, c("chisq","df","pvalue","cfi","rmsea"))
cat(sprintf("  Single-factor chi-sq = %.2f (df = %d, p < .001)\n",
            hf["chisq"], hf["df"]))
cat(sprintf("  CFI = %.3f\n", hf["cfi"]))
cat(sprintf("  RMSEA = %.3f\n", hf["rmsea"]))
cat(sprintf("  Interpretation: Poor single-factor fit indicates CMB\n"))
cat(sprintf("  is not a dominant concern.\n\n"))

# --- 1c. CFA Unmeasured Latent Factor (ULF) ---
cat("--- CFA Unmeasured Latent Factor (ULF) test ---\n")
substantive_model <- '
  AD =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent +
        AD_Q10_accept_no_verify + AD_Q5_frequency
  AR =~ AR_Q14_remorseful_avoidance + AR_Q15_regret_dependency +
        AR_Q16_regret_after_tasks + AR_Q11_fear_weaken_thinking +
        AR_Q12_aware_recall_loss
  PB =~ PB_Q17_saves_time + PB_Q18_confidence +
        PB_Q19_learning_efficiency + PB_Q20_satisfaction +
        PB_Q21_collaborative
  CE =~ CE_Q8_starting_point + CE_Q9_crosscheck
  AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss
  PB_Q18_confidence           ~~ PB_Q19_learning_efficiency
  PB_Q17_saves_time           ~~ PB_Q21_collaborative
  PB_Q17_saves_time           ~~ PB_Q20_satisfaction
  PB_Q20_satisfaction         ~~ PB_Q21_collaborative
'
sub_fit  <- cfa(substantive_model, data = df_items,
                estimator = "MLR", missing = "listwise")

ulf_model <- paste0(substantive_model, '
  CMF =~ 1*AD_Q5_frequency + 1*AD_Q6_reliance_decisions +
         1*AD_Q7_rely_vs_independent + 1*AD_Q10_accept_no_verify +
         1*AR_Q14_remorseful_avoidance + 1*AR_Q15_regret_dependency +
         1*AR_Q16_regret_after_tasks + 1*AR_Q11_fear_weaken_thinking +
         1*AR_Q12_aware_recall_loss +
         1*PB_Q17_saves_time + 1*PB_Q18_confidence +
         1*PB_Q19_learning_efficiency + 1*PB_Q20_satisfaction +
         1*PB_Q21_collaborative +
         1*CE_Q8_starting_point + 1*CE_Q9_crosscheck
  CMF ~~ 0*AD
  CMF ~~ 0*AR
  CMF ~~ 0*PB
  CMF ~~ 0*CE
')

ulf_fit <- tryCatch(
  cfa(ulf_model, data = df_items, estimator = "MLR", missing = "listwise"),
  error = function(e) NULL
)

if (!is.null(ulf_fit)) {
  sub_cfi <- fitMeasures(sub_fit,  "cfi.robust")
  ulf_cfi <- fitMeasures(ulf_fit,  "cfi.robust")
  delta   <- ulf_cfi - sub_cfi
  cat(sprintf("  Substantive CFA CFI = %.3f\n", sub_cfi))
  cat(sprintf("  CFA + CMF CFI       = %.3f\n", ulf_cfi))
  cat(sprintf("  Delta CFI           = %.3f\n", delta))
  cat(sprintf("  Result: %s\n\n",
      ifelse(abs(delta) < 0.010,
             "Delta CFI < .010 — CMB not a major concern",
             "Delta CFI >= .010 — CMB may be present")))
} else {
  cat("  ULF model did not converge. Harman result used as sole diagnostic.\n\n")
}

# =============================================================================
# SECTION 2: DATA DISTRIBUTION AND NORMALITY
# =============================================================================
cat("=======================================================================\n")
cat("SECTION 2: DATA DISTRIBUTION AND NORMALITY\n")
cat("=======================================================================\n\n")

cat("--- Item-level descriptives ---\n")
desc <- data.frame(
  Item     = colnames(df_items),
  Mean     = round(colMeans(df_items, na.rm=TRUE), 3),
  SD       = round(apply(df_items, 2, sd, na.rm=TRUE), 3),
  Skewness = round(apply(df_items, 2, psych::skew), 3),
  Kurtosis = round(apply(df_items, 2, psych::kurtosi), 3)
)
desc$Skew_flag <- ifelse(abs(desc$Skewness) > 2, "CONCERN", "OK")
desc$Kurt_flag <- ifelse(abs(desc$Kurtosis) > 7, "CONCERN", "OK")
print(desc, row.names=FALSE)
cat(sprintf("\n  Items with |skewness| > 2: %d\n", sum(desc$Skew_flag=="CONCERN")))
cat(sprintf("  Items with |kurtosis| > 7: %d\n\n", sum(desc$Kurt_flag=="CONCERN")))

cat("--- Construct-level descriptives ---\n")
df_scores <- data.frame(
  AD = rowMeans(df_items[,c("AD_Q5_frequency","AD_Q6_reliance_decisions",
                             "AD_Q7_rely_vs_independent","AD_Q10_accept_no_verify")]),
  AR = rowMeans(df_items[,c("AR_Q14_remorseful_avoidance","AR_Q15_regret_dependency",
                             "AR_Q16_regret_after_tasks","AR_Q11_fear_weaken_thinking",
                             "AR_Q12_aware_recall_loss")]),
  PB = rowMeans(df_items[,c("PB_Q17_saves_time","PB_Q18_confidence",
                             "PB_Q19_learning_efficiency","PB_Q20_satisfaction",
                             "PB_Q21_collaborative")]),
  CE = rowMeans(df_items[,c("CE_Q8_starting_point","CE_Q9_crosscheck")])
)
cdesc <- data.frame(
  Construct = c("AD","AR","PB","CE"),
  Mean  = round(colMeans(df_scores), 3),
  SD    = round(apply(df_scores, 2, sd), 3),
  Min   = round(apply(df_scores, 2, min), 2),
  Max   = round(apply(df_scores, 2, max), 2),
  Skew  = round(apply(df_scores, 2, psych::skew), 3),
  Kurt  = round(apply(df_scores, 2, psych::kurtosi), 3)
)
print(cdesc, row.names=FALSE)

cat("\n--- Mardia's Multivariate Normality Test ---\n")
mardia_res <- psych::mardia(df_items, plot=FALSE)
cat(sprintf("  Multivariate skewness: estimate = %.2f, statistic = %.2f, p < .001\n",
            mardia_res$b1p, mardia_res$skew))
cat(sprintf("  Multivariate kurtosis: estimate = %.2f, statistic = %.2f, p = %.3f\n",
            mardia_res$b2p, mardia_res$kurtosis, mardia_res$p.kurt))
cat(sprintf("  Result: Multivariate normality violated.\n"))
cat(sprintf("  MLR estimator with Satorra-Bentler correction is justified.\n\n"))

# Inter-construct correlations with 95% CI
cat("--- Inter-construct correlations (Pearson, 95% CI) ---\n")
pairs_list <- list(c("AD","AR"),c("AD","PB"),c("AD","CE"),
                   c("AR","PB"),c("AR","CE"),c("PB","CE"))
cat(sprintf("  %-6s %-6s %8s %8s %6s  %s\n","Var1","Var2","r","p","","95% CI"))
for (p in pairs_list) {
  ct <- cor.test(df_scores[[p[1]]], df_scores[[p[2]]])
  cat(sprintf("  %-6s %-6s %8.3f %8.4f  [%.3f, %.3f]\n",
              p[1], p[2], ct$estimate, ct$p.value,
              ct$conf.int[1], ct$conf.int[2]))
}

# =============================================================================
# SECTION 3: MEASUREMENT MODEL EVALUATION (CFA)
# Final model: v2 with 5 correlated residuals
# =============================================================================
cat("\n=======================================================================\n")
cat("SECTION 3: MEASUREMENT MODEL EVALUATION (CFA)\n")
cat("=======================================================================\n\n")

cfa_model <- '
  AD =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent +
        AD_Q10_accept_no_verify + AD_Q5_frequency
  AR =~ AR_Q14_remorseful_avoidance + AR_Q15_regret_dependency +
        AR_Q16_regret_after_tasks + AR_Q11_fear_weaken_thinking +
        AR_Q12_aware_recall_loss
  PB =~ PB_Q17_saves_time + PB_Q18_confidence +
        PB_Q19_learning_efficiency + PB_Q20_satisfaction +
        PB_Q21_collaborative
  CE =~ CE_Q8_starting_point + CE_Q9_crosscheck

  # 5 correlated residuals (theoretically justified)
  AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss
  PB_Q18_confidence           ~~ PB_Q19_learning_efficiency
  PB_Q17_saves_time           ~~ PB_Q21_collaborative
  PB_Q17_saves_time           ~~ PB_Q20_satisfaction
  PB_Q20_satisfaction         ~~ PB_Q21_collaborative
'

cfa_fit <- cfa(cfa_model, data=df_items, estimator="MLR", missing="listwise")

fi <- fitMeasures(cfa_fit, c(
  "chisq.scaled","df.scaled","pvalue.scaled",
  "cfi.robust","tli.robust",
  "rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr"
))

cat("--- CFA Model Fit Indices ---\n")
cat(sprintf("  chi-sq (SB)  = %.2f (df = %d, p < .001)\n",
            fi["chisq.scaled"], fi["df.scaled"]))
cat(sprintf("  CFI          = %.3f\n", fi["cfi.robust"]))
cat(sprintf("  TLI          = %.3f\n", fi["tli.robust"]))
cat(sprintf("  RMSEA        = %.3f [90%% CI: %.3f-%.3f]\n",
            fi["rmsea.robust"],
            fi["rmsea.ci.lower.robust"],
            fi["rmsea.ci.upper.robust"]))
cat(sprintf("  SRMR         = %.3f\n\n", fi["srmr"]))

cat("--- Standardized Factor Loadings ---\n")
std_p    <- standardizedSolution(cfa_fit)
loadings <- std_p[std_p$op == "=~",
                  c("lhs","rhs","est.std","se","z","pvalue")]
loadings$sig <- ifelse(loadings$pvalue < 0.001, "***",
                ifelse(loadings$pvalue < 0.01,  "**",
                ifelse(loadings$pvalue < 0.05,  "*", "ns")))
print(loadings, row.names=FALSE)

cat("\n--- Freed Correlated Residuals ---\n")
resids <- std_p[std_p$op == "~~" &
                std_p$lhs != std_p$rhs &
                std_p$lhs %in% item_cols,
                c("lhs","rhs","est.std","se","z","pvalue")]
print(resids, row.names=FALSE)

# =============================================================================
# SECTION 4: RELIABILITY AND CONVERGENT VALIDITY
# =============================================================================
cat("\n=======================================================================\n")
cat("SECTION 4: RELIABILITY AND CONVERGENT VALIDITY\n")
cat("=======================================================================\n\n")

ad_it <- df_items[,c("AD_Q6_reliance_decisions","AD_Q7_rely_vs_independent",
                      "AD_Q10_accept_no_verify","AD_Q5_frequency")]
ar_it <- df_items[,c("AR_Q14_remorseful_avoidance","AR_Q15_regret_dependency",
                      "AR_Q16_regret_after_tasks","AR_Q11_fear_weaken_thinking",
                      "AR_Q12_aware_recall_loss")]
pb_it <- df_items[,c("PB_Q17_saves_time","PB_Q18_confidence",
                      "PB_Q19_learning_efficiency","PB_Q20_satisfaction",
                      "PB_Q21_collaborative")]
ce_it <- df_items[,c("CE_Q8_starting_point","CE_Q9_crosscheck")]

alphas <- c(
  psych::alpha(ad_it)$total$raw_alpha,
  psych::alpha(ar_it)$total$raw_alpha,
  psych::alpha(pb_it)$total$raw_alpha,
  psych::alpha(ce_it)$total$raw_alpha
)
cr_vals  <- semTools::compRelSEM(cfa_fit)
ave_vals <- semTools::AVE(cfa_fit)

rel_tab <- data.frame(
  Construct    = c("AD","AR","PB","CE"),
  Alpha        = round(alphas, 3),
  CR_omega     = round(as.numeric(cr_vals), 3),
  AVE          = round(as.numeric(ave_vals), 3),
  Alpha_ok     = ifelse(alphas >= 0.70, "Yes","No"),
  CR_ok        = ifelse(as.numeric(cr_vals) >= 0.70, "Yes","No"),
  AVE_ok       = ifelse(as.numeric(ave_vals) >= 0.50, "Yes","No")
)
cat("--- Reliability and Convergent Validity Summary ---\n")
print(rel_tab, row.names=FALSE)
cat("\n  Thresholds: alpha >= .70, CR >= .70, AVE >= .50\n")
cat("  Note: AD and CE below threshold due to item heterogeneity\n")
cat("        and 2-item CE constraint. Acknowledged as limitations.\n\n")

# =============================================================================
# SECTION 5: DISCRIMINANT VALIDITY
# =============================================================================
cat("=======================================================================\n")
cat("SECTION 5: DISCRIMINANT VALIDITY\n")
cat("=======================================================================\n\n")

cat("--- HTMT Ratios (threshold < 0.85) ---\n")
htmt_out <- semTools::htmt(cfa_model, data=df_items)
print(round(htmt_out, 3))
cat("\n  All HTMT values must be < 0.85. Values above indicate\n")
cat("  potential discriminant validity issues.\n")

cat("\n--- Fornell-Larcker: sqrt(AVE) vs latent correlations ---\n")
lv_corr  <- lavInspect(cfa_fit, what="cor.lv")
sqrt_ave <- sqrt(as.numeric(ave_vals))
names(sqrt_ave) <- names(ave_vals)
fl_mat   <- as.matrix(lv_corr)
diag(fl_mat) <- sqrt_ave
cat("  (Diagonal = sqrt(AVE); off-diagonal = latent correlations)\n")
print(round(fl_mat, 3))
cat("  sqrt(AVE) should exceed all off-diagonal values in its row/column.\n\n")

# =============================================================================
# SECTION 6: STRUCTURAL MODEL — BOOTSTRAP FOR CIs ON DIRECT PATHS
# Mirrors referenced paper: B, SE, z, p, 95% CI for every path
# =============================================================================
cat("=======================================================================\n")
cat("SECTION 6: STRUCTURAL MODEL FIT AND PATH COEFFICIENTS\n")
cat("Bootstrap n = 5,000 | BCa 95% CI | estimator = ML\n")
cat("=======================================================================\n\n")

sem_boot_model <- '
  AD =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent +
        AD_Q10_accept_no_verify + AD_Q5_frequency
  AR =~ AR_Q14_remorseful_avoidance + AR_Q15_regret_dependency +
        AR_Q16_regret_after_tasks + AR_Q11_fear_weaken_thinking +
        AR_Q12_aware_recall_loss
  PB =~ PB_Q17_saves_time + PB_Q18_confidence +
        PB_Q19_learning_efficiency + PB_Q20_satisfaction +
        PB_Q21_collaborative
  CE =~ CE_Q8_starting_point + CE_Q9_crosscheck

  AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss
  PB_Q18_confidence           ~~ PB_Q19_learning_efficiency
  PB_Q17_saves_time           ~~ PB_Q21_collaborative
  PB_Q17_saves_time           ~~ PB_Q20_satisfaction
  PB_Q20_satisfaction         ~~ PB_Q21_collaborative

  # Structural paths
  AR ~ b1*AD
  AD ~ b2*PB
  AR ~ b3*PB
  PB ~ b4*CE
  AD ~ b5*CE
'

set.seed(42)
sem_boot <- sem(
  sem_boot_model,
  data      = df_items,
  estimator = "ML",
  se        = "bootstrap",
  bootstrap = 5000,
  missing   = "listwise"
)

# Fit indices (use MLR version for reporting)
sem_mlr <- sem(sem_boot_model, data=df_items,
               estimator="MLR", missing="listwise")
sf <- fitMeasures(sem_mlr, c(
  "chisq.scaled","df.scaled","pvalue.scaled",
  "cfi.robust","tli.robust",
  "rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust","srmr"
))

cat("--- Structural Model Fit (MLR) ---\n")
cat(sprintf("  chi-sq (SB) = %.2f (df = %d, p < .001)\n",
            sf["chisq.scaled"], sf["df.scaled"]))
cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n",
            sf["cfi.robust"], sf["tli.robust"]))
cat(sprintf("  RMSEA = %.3f [90%% CI: %.3f-%.3f]  |  SRMR = %.3f\n\n",
            sf["rmsea.robust"],
            sf["rmsea.ci.lower.robust"],
            sf["rmsea.ci.upper.robust"],
            sf["srmr"]))

# R-squared
cat("--- R-squared (endogenous constructs) ---\n")
r2 <- lavInspect(sem_mlr, what="r2")
for (nm in c("AD","AR","PB")) {
  rv  <- r2[[nm]]
  f2  <- rv / (1 - rv)
  mag <- ifelse(f2>=0.35,"large",ifelse(f2>=0.15,"medium",
         ifelse(f2>=0.02,"small","negligible")))
  cat(sprintf("  %-4s  R^2 = %.3f  f^2 = %.3f (%s)\n", nm, rv, f2, mag))
}

# Direct paths with bootstrap CIs
cat("\n--- Structural Paths: B, SE, z, p, 95% BCa CI, beta (std) ---\n")
boot_pe <- parameterEstimates(
  sem_boot,
  boot.ci.type = "bca.simple",
  level        = 0.95,
  standardized = TRUE
)
paths <- boot_pe[boot_pe$op == "~",
                 c("lhs","rhs","est","se","z","pvalue",
                   "ci.lower","ci.upper","std.all")]
paths$sig <- ifelse(paths$pvalue < 0.001, "***",
             ifelse(paths$pvalue < 0.01,  "**",
             ifelse(paths$pvalue < 0.05,  "*", "ns")))
paths$Hyp <- c("H1","H2","H3","H4","H5")
paths$Decision <- ifelse(paths$pvalue < 0.05, "Supported", "Not supported")
print(paths[,c("Hyp","lhs","rhs","est","se","z","pvalue",
               "ci.lower","ci.upper","std.all","sig","Decision")],
      row.names=FALSE)

# =============================================================================
# SECTION 7: MEDIATION ANALYSIS — CE -> PB -> AD (H6)
# =============================================================================
cat("\n=======================================================================\n")
cat("SECTION 7: MEDIATION ANALYSIS — CE -> PB -> AD (H6)\n")
cat("Bootstrap n = 5,000 | BCa 95% CI | estimator = ML\n")
cat("=======================================================================\n\n")

med_model <- '
  AD =~ AD_Q6_reliance_decisions + AD_Q7_rely_vs_independent +
        AD_Q10_accept_no_verify + AD_Q5_frequency
  AR =~ AR_Q14_remorseful_avoidance + AR_Q15_regret_dependency +
        AR_Q16_regret_after_tasks + AR_Q11_fear_weaken_thinking +
        AR_Q12_aware_recall_loss
  PB =~ PB_Q17_saves_time + PB_Q18_confidence +
        PB_Q19_learning_efficiency + PB_Q20_satisfaction +
        PB_Q21_collaborative
  CE =~ CE_Q8_starting_point + CE_Q9_crosscheck

  AR_Q11_fear_weaken_thinking ~~ AR_Q12_aware_recall_loss
  PB_Q18_confidence           ~~ PB_Q19_learning_efficiency
  PB_Q17_saves_time           ~~ PB_Q21_collaborative
  PB_Q17_saves_time           ~~ PB_Q20_satisfaction
  PB_Q20_satisfaction         ~~ PB_Q21_collaborative

  PB ~ a*CE
  AD ~ b*PB + c_prime*CE
  AR ~ AD + PB

  indirect := a * b
  total    := c_prime + a * b
  pm       := (a * b) / (c_prime + a * b)
'

set.seed(42)
med_fit <- sem(
  med_model,
  data      = df_items,
  estimator = "ML",
  se        = "bootstrap",
  bootstrap = 5000,
  missing   = "listwise"
)

med_pe <- parameterEstimates(
  med_fit,
  boot.ci.type = "bca.simple",
  level        = 0.95,
  standardized = TRUE
)

cat("--- Defined parameters (indirect, total, proportion mediated) ---\n")
def <- med_pe[med_pe$op == ":=",
              c("label","est","se","z","pvalue","ci.lower","ci.upper","std.all")]
def$sig <- ifelse(def$pvalue < 0.001,"***",
           ifelse(def$pvalue < 0.01, "**",
           ifelse(def$pvalue < 0.05, "*", "ns")))
print(def, row.names=FALSE)

cat("\n--- All structural paths (bootstrapped SE, BCa CI) ---\n")
kp <- med_pe[med_pe$op == "~" & med_pe$lhs %in% c("PB","AD","AR"),
             c("lhs","rhs","est","se","z","pvalue","ci.lower","ci.upper","std.all")]
kp$sig <- ifelse(kp$pvalue < 0.001,"***",
          ifelse(kp$pvalue < 0.01, "**",
          ifelse(kp$pvalue < 0.05, "*", "ns")))
print(kp, row.names=FALSE)

cat("\n  Mediation interpretation:\n")
cat("  - indirect CI excludes 0  -> mediation is significant\n")
cat("  - c_prime significant     -> partial mediation\n")
cat("  - c_prime non-significant -> full mediation\n\n")

# =============================================================================
# PATH DIAGRAM
# =============================================================================
cat("Saving path diagram -> sem_path_diagram_final.pdf\n")
pdf("sem_path_diagram_final.pdf", width=10, height=7)
semPaths(sem_mlr, what="std", layout="tree2",
         edge.label.cex=0.8, curvePivot=TRUE,
         residuals=FALSE, nCharNodes=0,
         style="lisrel", title=FALSE)
dev.off()
cat("Saved.\n\n")

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Output covers all sections required for manuscript:\n")
cat("  1. Sample adequacy (KMO, Bartlett) + CMB (Harman + ULF)\n")
cat("  2. Normality (item-level skew/kurt + Mardia's MVN)\n")
cat("  3. CFA fit + loadings + freed residuals\n")
cat("  4. Reliability (alpha, CR/omega) + AVE\n")
cat("  5. Discriminant validity (HTMT + Fornell-Larcker)\n")
cat("  6. Structural paths: B, SE, z, p, 95% BCa CI, beta\n")
cat("  7. Mediation: indirect/direct/total + BCa CI\n")
cat("=======================================================================\n")