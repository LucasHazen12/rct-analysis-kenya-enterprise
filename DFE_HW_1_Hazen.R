# DFE HW 1
# Lucas Hazen

setwd("/Users/lucashazen/Desktop/R Materials")

# install.packages(c("tidyverse", "WebPower", "ICC", "fishmethods",
#                    "stargazer", "modelsummary", "AER", "sandwich", "lmtest",
#                    "margins", "kableExtra", "MatchIt"))

# Load libraries
library(tidyverse)
library(WebPower)
library(ICC)
library(fishmethods)
library(stargazer)
library(modelsummary)
library(margins)
library(kableExtra)
library(MatchIt)

# Load Data

df <- read.csv("HW1.csv")

# Convert 'trusting' from "yes"/"no" string to a binary 0/1
df <- df %>%
  mutate(trusting = ifelse(trusting == "yes", 1, 0))

# Q1: Power Curve for Cluster-Randomized Trial
# Use round 1 (baseline) data only

r1 <- df %>% filter(round == 1)

# Step 1: Calculate J (number of clusters)
J <- n_distinct(r1$vill_id)
cat("Number of clusters (J):", J, "\n")

# Step 2: Calculate n (average number of units per cluster)
n_per_cluster <- r1 %>%
  group_by(vill_id) %>%
  summarise(n = n())

n_avg <- mean(n_per_cluster$n)
cat("Average units per cluster (n):", round(n_avg, 2), "\n")

# Step 3: Calculate ICC using the ICC package
icc_result <- ICCbare(x = factor(vill_id), y = earnings, data = r1)
cat("Intra-cluster correlation (ICC):", round(icc_result, 4), "\n")

# Step 4: R-squared of earnings on baseline covariates
# This the covariate-adjusted scenario
r1_cov <- r1 %>% drop_na(earnings, trusting, risk_amount, age, urban, highest_grade)

cov_model <- lm(earnings ~ trusting + risk_amount + age + urban + highest_grade,
                data = r1_cov)
r2_covariates <- summary(cov_model)$r.squared
cat("R-squared from baseline covariates:", round(r2_covariates, 4), "\n")

# Step 5: SD of earnings at baseline (needed for MDE in Ksh)
sd_earnings <- sd(r1$earnings, na.rm = TRUE)
cat("SD of baseline earnings (Ksh):", round(sd_earnings, 2), "\n")

# Step 6: Plot Power Curve for CRT
# Two scenarios: no covariates (R2=0) vs. with covariates (R2=0.0879)

# Effect sizes to loop over (in SD units)
effect_sizes <- seq(0.01, 1.0, by = 0.01)

# Storage
power_no_cov  <- numeric(length(effect_sizes))
power_with_cov <- numeric(length(effect_sizes))

for (i in seq_along(effect_sizes)) {
  f <- effect_sizes[i]
  
  # Scenario 1: No covariates (R2 = 0)
  power_no_cov[i] <- wp.crt2arm(
    f     = f,
    n     = n_avg,
    J     = J,
    icc   = icc_result,
    alpha = 0.05,
    power = NULL
  )$power
  
  # Scenario 2: With covariates (R2 = 0.0879)
  # Covariates reduce residual variance: effective f is scaled up
  f_adj <- f / sqrt(1 - r2_covariates)
  power_with_cov[i] <- wp.crt2arm(
    f     = f_adj,
    n     = n_avg,
    J     = J,
    icc   = icc_result,
    alpha = 0.05,
    power = NULL
  )$power
}

# Build data frame for plotting
power_df <- data.frame(
  effect_size   = rep(effect_sizes, 2),
  power         = c(power_no_cov, power_with_cov),
  scenario      = rep(c("No Covariates (R²=0)", "With Covariates (R²=0.088)"), 
                      each = length(effect_sizes))
)

# Plot
ggplot(power_df, aes(x = effect_size, y = power, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  annotate("text", x = 0.85, y = 0.82, label = "Power = 0.80", 
           color = "red", size = 3.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Power Curve: Cluster-Randomized Trial (Kenya Enterprise Program)",
    subtitle = paste0("J=168 villages, n≈4.79 per cluster, ICC=0.056, α=0.05"),
    x        = "Effect Size (in SD units)",
    y        = "Statistical Power",
    linetype = "Scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("power_curve_HW1.png", width = 7, height = 5, units = "in")


# Q2: Minimum Detectable Effect (MDE) at 80% Power
# Solve for f (effect size in SD units) where power = 0.80

# Scenario 1: No covariates (R2 = 0)
mde_no_cov <- wp.crt2arm(
  f     = NULL,
  n     = n_avg,
  J     = J,
  icc   = icc_result,
  alpha = 0.05,
  power = 0.8
)$f

# Scenario 2: With covariates — covariates reduce residual SD by sqrt(1 - R2)
# So the MDE in original SD units is scaled back down
mde_with_cov_adj <- wp.crt2arm(
  f     = NULL,
  n     = n_avg,
  J     = J,
  icc   = icc_result,
  alpha = 0.05,
  power = 0.8
)$f * sqrt(1 - r2_covariates)

# Convert from SD units to Kenyan Shillings
mde_no_cov_ksh   <- mde_no_cov * sd_earnings
mde_with_cov_ksh <- mde_with_cov_adj * sd_earnings

cat("=== MDE Results ===\n")
cat("No covariates:   ", round(mde_no_cov, 4), "SD units =", 
    round(mde_no_cov_ksh, 0), "Ksh\n")
cat("With covariates: ", round(mde_with_cov_adj, 4), "SD units =", 
    round(mde_with_cov_ksh, 0), "Ksh\n")
cat("Difference in MDE (Ksh):", 
    round(mde_no_cov_ksh - mde_with_cov_ksh, 0), "Ksh\n")


# Q3: Features of the Treatment

# Q3a: Unit of intervention
village_offer_r1 <- df %>%
  filter(round == 1) %>%
  group_by(vill_id) %>%
  summarise(n_offered_vals = n_distinct(offered))

cat("Villages with mixed 'offered' values (round 1 only):", 
    sum(village_offer_r1$n_offered_vals > 1), "\n")

# Q3b: Balance table on baseline covariates only
# Offered is recorded in round 2 only, merge r1 covariates with r2 treatment status
r1_covs <- df %>%
  filter(round == 1) %>%
  dplyr::select(indiv_id, earnings, trusting, risk_amount, age, urban, highest_grade)

r2_treatment <- df %>%
  filter(round == 2) %>%
  dplyr::select(indiv_id, offered, treated)

r1_balance <- r1_covs %>%
  left_join(r2_treatment, by = "indiv_id") %>%
  filter(!is.na(offered)) %>%
  mutate(offered = as.factor(offered))

covs <- c("earnings", "trusting", "risk_amount", "age", "urban", "highest_grade")

balance_results <- lapply(covs, function(var) {
  t_result <- t.test(as.formula(paste(var, "~ offered")), data = r1_balance)
  data.frame(
    Variable     = var,
    Mean_Control = round(t_result$estimate[1], 3),
    Mean_Treated = round(t_result$estimate[2], 3),
    Difference   = round(t_result$estimate[2] - t_result$estimate[1], 3),
    P_value      = round(t_result$p.value, 3)
  )
})

balance_table <- do.call(rbind, balance_results)
rownames(balance_table) <- NULL
print(balance_table)

# Export balance table as HTML using stargazer
stargazer(as.data.frame(balance_table),
          type = "html",
          out = "Q3b_balance_table.html",
          summary = FALSE,
          rownames = FALSE,
          title = "Balance of Randomization on Baseline Covariates")

cat("Balance table saved to Q3b_balance_table.html\n")

# Q3c: One-sided or two-sided non-compliance?
r2 <- df %>% filter(round == 2)

control_treated  <- r2 %>% filter(offered == 0, treated == 1) %>% nrow()
treat_not_treated <- r2 %>% filter(offered == 1, treated == 0) %>% nrow()

cat("\nControl units who received treatment:", control_treated, "\n")
cat("Offered units who did NOT take treatment:", treat_not_treated, "\n")

# Q3d: Uptake rate
uptake_rate <- r2 %>%
  filter(offered == 1) %>%
  summarise(uptake = mean(treated, na.rm = TRUE))

cat("Uptake rate among offered:", round(uptake_rate$uptake, 4),
    "(", round(uptake_rate$uptake * 100, 1), "%)\n")


# Q4: Probit Regression on Program Uptake

# Build probit dataset
# Treatment outcome from round 2, covariates from round 1
r1_probit <- df %>%
  filter(round == 1) %>%
  dplyr::select(indiv_id, trusting, risk_amount, age, urban, highest_grade)

r2_probit <- df %>%
  filter(round == 2, offered == 1) %>%
  dplyr::select(indiv_id, treated)

probit_data <- r2_probit %>%
  left_join(r1_probit, by = "indiv_id")

cat("Probit sample size (offered=1):", nrow(probit_data), "\n")
cat("Treated vs not-treated in sample:\n")
print(table(probit_data$treated))

# Run probit
probit_model <- glm(
  treated ~ trusting + risk_amount + age + urban + highest_grade,
  data   = probit_data,
  family = binomial(link = "probit")
)

summary(probit_model)

# Marginal effects (easier to interpret than raw coefficients)
marginal_effects <- margins(probit_model)
summary(marginal_effects)


# Q5: Intention to Treat (ITT) Effect

# ITT uses 'offered' as the treatment indicator (regardless of takeup)
# Use round 2 (endline) earnings as the outcome
# Covariates are time-invariant, measured at baseline (round 1)

# Build the ITT dataset: endline outcome + treatment + baseline covariates
r1_itt <- df %>%
  filter(round == 1) %>%
  dplyr::select(indiv_id, earnings_baseline = earnings,
                trusting, risk_amount, age, urban, highest_grade)

r2_itt <- df %>%
  filter(round == 2) %>%
  dplyr::select(indiv_id, earnings_endline = earnings, offered)

itt_data <- r2_itt %>%
  left_join(r1_itt, by = "indiv_id")

cat("ITT sample size:", nrow(itt_data), "\n")

# --- Q5a: ITT with no covariates ---
itt_simple <- lm(earnings_endline ~ offered, data = itt_data)
summary(itt_simple)

# --- Q5b: ANCOVA — endline on offered + baseline outcome + covariates ---
itt_ancova <- lm(earnings_endline ~ offered + earnings_baseline +
                   trusting + risk_amount + age + urban + highest_grade,
                 data = itt_data)
summary(itt_ancova)

# --- Display both side by side ---
stargazer(itt_simple, itt_ancova,
          type = "text",
          title = "ITT Estimates: Simple OLS vs ANCOVA",
          column.labels = c("No Covariates", "ANCOVA"),
          covariate.labels = c("Offered Treatment", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Export stargazer table as HTML
stargazer(itt_simple, itt_ancova,
          type = "html",
          out = "Q5_ITT_table.html",
          title = "ITT Estimates: Simple OLS vs ANCOVA",
          column.labels = c("No Covariates", "ANCOVA"),
          covariate.labels = c("Offered Treatment", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Q6: Randomization Inference (500 repetitions)
# Mimicking village-level assignment (the true assignment rule)

library(sandwich)
library(lmtest)

set.seed(123)
n_sims <- 500

# Get village-level treatment status from round 2
village_treatment <- df %>%
  filter(round == 2) %>%
  group_by(vill_id) %>%
  summarise(offered = first(offered)) %>%
  ungroup()

n_treated_villages <- sum(village_treatment$offered)
n_total_villages   <- nrow(village_treatment)

cat("Total villages:", n_total_villages, "\n")
cat("Treated villages:", n_treated_villages, "\n")

# Store simulated ITT coefficients
sim_itt <- numeric(n_sims)

for (i in 1:n_sims) {
  village_treatment$offered_sim <- 0
  village_treatment$offered_sim[
    sample(n_total_villages, n_treated_villages)] <- 1
  
  sim_data <- df %>%
    filter(round == 2) %>%
    left_join(village_treatment %>%
                dplyr::select(vill_id, offered_sim), by = "vill_id")
  
  sim_model <- lm(earnings ~ offered_sim, data = sim_data)
  sim_itt[i] <- coef(sim_model)["offered_sim"]
}

# RI standard error
ri_se <- sd(sim_itt)

# Rebuild itt_data with vill_id included for clustering
r1_itt <- df %>%
  filter(round == 1) %>%
  dplyr::select(indiv_id, earnings_baseline = earnings,
                trusting, risk_amount, age, urban, highest_grade)

r2_itt <- df %>%
  filter(round == 2) %>%
  dplyr::select(indiv_id, earnings_endline = earnings, offered)

itt_data <- r2_itt %>%
  left_join(r1_itt, by = "indiv_id") %>%
  left_join(df %>% filter(round == 2) %>%
              dplyr::select(indiv_id, vill_id), by = "indiv_id")

# Re-run itt_simple on updated itt_data
itt_simple <- lm(earnings_endline ~ offered, data = itt_data)

# Clustered SE
itt_clustered <- coeftest(itt_simple,
                          vcov = vcovCL(itt_simple,
                                        cluster = ~vill_id,
                                        data = itt_data))

parametric_se_plain     <- summary(itt_simple)$coefficients["offered", "Std. Error"]
parametric_se_clustered <- itt_clustered["offered", "Std. Error"]

cat("\nRandomization Inference SE:       ", round(ri_se, 2), "\n")
cat("Parametric SE (unclustered):      ", round(parametric_se_plain, 2), "\n")
cat("Parametric SE (clustered):        ", round(parametric_se_clustered, 2), "\n")

# Plot RI distribution
ri_df <- data.frame(sim_itt = sim_itt)
ggplot(ri_df, aes(x = sim_itt)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = coef(itt_simple)["offered"],
             color = "red", linewidth = 1, linetype = "dashed") +
  annotate("text", x = coef(itt_simple)["offered"] + 300,
           y = 30, label = "Observed ITT", color = "red", size = 3.5) +
  labs(title = "Randomization Inference Distribution (500 simulations)",
       subtitle = "Village-level randomization, no covariates",
       x = "Simulated ITT Coefficient", y = "Count") +
  theme_bw()

ggsave("RI_distribution_HW1.png", width = 7, height = 4, units = "in")


# Q7: IV Treatment on the Treated (ToT) Regression

# With one-sided non-compliance, OLS on 'treated' is biased
# IV solution: instrument 'treated' with 'offered' (the randomized assignment)
# This gives the LATE, the effect on compliers only

library(AER)

# Build dataset: same as ITT but need 'treated' as the endogenous variable
r1_late <- df %>%
  filter(round == 1) %>%
  dplyr::select(indiv_id, earnings_baseline = earnings,
                trusting, risk_amount, age, urban, highest_grade)

r2_late <- df %>%
  filter(round == 2) %>%
  dplyr::select(indiv_id, earnings_endline = earnings, offered, treated)

late_data <- r2_late %>%
  left_join(r1_late, by = "indiv_id")

# Q7a: Simple IV (no covariates)
iv_simple <- ivreg(earnings_endline ~ treated | offered,
                   data = late_data)
summary(iv_simple)

# Q7b: IV with ANCOVA controls
iv_ancova <- ivreg(earnings_endline ~ treated + earnings_baseline +
                     trusting + risk_amount + age + urban + highest_grade |
                     offered + earnings_baseline +
                     trusting + risk_amount + age + urban + highest_grade,
                   data = late_data)
summary(iv_ancova)

# Display both side by side 
stargazer(iv_simple, iv_ancova,
          type = "text",
          title = "ToT Estimates: IV (2SLS)",
          column.labels = c("No Covariates", "With Covariates"),
          covariate.labels = c("Treated", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Sanity check: LATE = ITT / Uptake rate
itt_estimate <- coef(itt_simple)["offered"]
uptake       <- uptake_rate$uptake
late_manual  <- itt_estimate / uptake

cat("\nManual LATE check (ITT / uptake rate):", round(late_manual, 0), "Ksh\n")
cat("IV LATE estimate:", round(coef(iv_simple)["treated"], 0), "Ksh\n")

# Q7: Export regression table as HTML for write-up
stargazer(iv_simple, iv_ancova,
          type = "html",
          out = "Q7_ToT_IV_table.html",
          title = "ToT Estimates: IV (2SLS)",
          column.labels = c("No Covariates", "With Covariates"),
          covariate.labels = c("Treated", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Q8: Matching Estimator of ToT

# Matching estimates the ToT by finding control units that look
# like treated units on observable characteristics, then comparing
# outcomes. Assumption: selection on observables (no unobserved
# confounders drive both treatment take-up and outcomes)

# Build matching dataset: round 2 outcomes, round 1 covariates
# Subsample: only those offered treatment (same as probit in Q4)
# Match treated vs non-treated among the offered group

match_data <- late_data %>%
  filter(!is.na(treated)) %>%
  mutate(treated = as.integer(treated))

cat("Matching sample size:", nrow(match_data), "\n")
cat("Treated:", sum(match_data$treated), "\n")
cat("Control (offered but not treated):", sum(match_data$treated == 0), "\n")

# Nearest neighbor matching on all baseline covariates
# Using 1-to-1 nearest neighbor matching with replacement
m_out <- matchit(treated ~ earnings_baseline + trusting + risk_amount +
                   age + urban + highest_grade,
                 data   = match_data,
                 method = "nearest",
                 ratio  = 1,
                 replace = TRUE)

summary(m_out)

# Extract matched data and estimate ToT
matched_data <- match.data(m_out)

tot_match <- lm(earnings_endline ~ treated + earnings_baseline +
                  trusting + risk_amount + age + urban + highest_grade,
                data    = matched_data,
                weights = weights)

summary(tot_match)

# Display results
stargazer(tot_match,
          type = "text",
          title = "ToT Estimate: Matching Estimator",
          covariate.labels = c("Treated", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Export as HTML
stargazer(tot_match,
          type = "html",
          out  = "Q8_matching_table.html",
          title = "ToT Estimate: Matching Estimator",
          covariate.labels = c("Treated", "Baseline Earnings",
                               "Trusting", "Risk Amount", "Age",
                               "Urban", "Highest Grade"),
          dep.var.labels = "Endline Earnings (Ksh)",
          omit.stat = c("f", "ser"),
          digits = 2)

# Compare IV vs Matching ToT
cat("\nIV ToT estimate:      ", round(coef(iv_simple)["treated"], 0), "Ksh\n")
cat("Matching ToT estimate:", round(coef(tot_match)["treated"], 0), "Ksh\n")

# Export balance plot for write-up
plot(summary(m_out), main = "Covariate Balance Before and After Matching")


# Q10: Heterogeneous Intention to Treat Effects

# Test whether the ITT effect varies across subgroups
# by interacting 'offered' with each baseline covariate
# Concern: multiple inference, testing many subgroups
# inflates the chance of false positives
# Use the itt_data which already has vill_id for clustering

# Median split each continuous variable for subgroup analysis
itt_data <- itt_data %>%
  mutate(
    above_median_earnings  = as.integer(earnings_baseline > median(earnings_baseline, na.rm = TRUE)),
    above_median_age       = as.integer(age > median(age, na.rm = TRUE)),
    above_median_risk      = as.integer(risk_amount > median(risk_amount, na.rm = TRUE)),
    above_median_grade     = as.integer(highest_grade > median(highest_grade, na.rm = TRUE))
  )

# Run ITT with interaction for each subgroup variable
# Binary variables: trusting, urban
# Median-split variables: earnings, age, risk_amount, highest_grade

het_vars <- list(
  "Trusting"         = "trusting",
  "Urban"            = "urban",
  "Above Median Baseline Earnings" = "above_median_earnings",
  "Above Median Age"               = "above_median_age",
  "Above Median Risk Amount"       = "above_median_risk",
  "Above Median Education"         = "above_median_grade"
)

het_results <- lapply(names(het_vars), function(label) {
  var <- het_vars[[label]]
  formula <- as.formula(paste("earnings_endline ~ offered *", var))
  model <- lm(formula, data = itt_data)
  
  # Get clustered SEs
  cl_se <- coeftest(model, vcov = vcovCL(model, cluster = ~vill_id, data = itt_data))
  
  # Interaction term is the heterogeneous effect
  int_term <- paste0("offered:", var)
  if (int_term %in% rownames(cl_se)) {
    data.frame(
      Subgroup    = label,
      Interaction = round(cl_se[int_term, "Estimate"], 1),
      SE          = round(cl_se[int_term, "Std. Error"], 1),
      P_value     = round(cl_se[int_term, "Pr(>|t|)"], 3)
    )
  }
})

het_table <- do.call(rbind, het_results)
rownames(het_table) <- NULL
print(het_table)

# Bonferroni correction for multiple inference
# Tested 6 hypotheses, so adjust p-values
het_table$P_bonferroni <- round(p.adjust(het_table$P_value, method = "bonferroni"), 3)
het_table$P_holm       <- round(p.adjust(het_table$P_value, method = "holm"), 3)

cat("\nHeterogeneous ITT results with multiple inference correction:\n")
print(het_table)

# Export as HTML
stargazer(as.data.frame(het_table),
          type = "html",
          out = "Q10_heterogeneous_ITT.html",
          summary = FALSE,
          rownames = FALSE,
          title = "Heterogeneous ITT Effects with Multiple Inference Correction")

cat("\nTable saved to Q10_heterogeneous_ITT.html\n")


# Q11: Lee Bounds and Attrition

# Q11a:
# If attrition is 0% in treatment and 20% in control:
# The treatment group is fully observed, control is trimmed
# Lee bounds will always contain the true effect because
# we only need to trim the larger group (control here)
# and the true effect is always within the bounds by construction

# Q11b: 20% attrition in treatment, 40% in control
# Artificially impose these attrition rates, then estimate
# the attrited ITT and Lee Bounds

# Use the itt_data (round 2 endline cross section)
# Treatment group: offered == 1, Control group: offered == 0

set.seed(42)

# Separate treatment and control
treat_group   <- itt_data %>% filter(offered == 1)
control_group <- itt_data %>% filter(offered == 0)

cat("Treatment group size:", nrow(treat_group), "\n")
cat("Control group size:", nrow(control_group), "\n")

# Artificially drop 20% of treatment and 40% of control (random attrition)
treat_attrited <- treat_group %>%
  slice_sample(prop = 0.80)   # keep 80% of treatment

control_attrited <- control_group %>%
  slice_sample(prop = 0.60)   # keep 60% of control

attrited_data <- bind_rows(treat_attrited, control_attrited)

cat("\nAfter attrition:\n")
cat("Treatment retained:", nrow(treat_attrited), "\n")
cat("Control retained:", nrow(control_attrited), "\n")

# Attrited ITT estimate
itt_attrited <- lm(earnings_endline ~ offered, data = attrited_data)
summary(itt_attrited)

attrited_coef <- coef(itt_attrited)["offered"]
cat("\nAttrited ITT estimate:", round(attrited_coef, 0), "Ksh\n")

# Lee Bounds
# Step 1: Calculate the trimming fraction
# Attrition is higher in control (40%) than treatment (20%)
# So trim from the treatment group
# Trimming fraction p = (s_t - s_c) / s_t
# where s_t = survival rate in treatment, s_c = survival rate in control

s_treat   <- 0.80   # survival rate in treatment
s_control <- 0.60   # survival rate in control

trim_fraction <- (s_treat - s_control) / s_treat
cat("\nTrimming fraction:", round(trim_fraction, 4), "\n")

# Step 2: Sort treatment group by outcome and trim
# Upper bound: drop bottom trim_fraction of treatment outcomes
# Lower bound: drop top trim_fraction of treatment outcomes

n_trim <- floor(trim_fraction * nrow(treat_attrited))
cat("Number of observations to trim:", n_trim, "\n")

treat_sorted <- treat_attrited %>%
  arrange(earnings_endline)

# Upper bound: drop lowest outcomes from treatment
treat_upper <- treat_sorted %>% slice((n_trim + 1):n())
# Lower bound: drop highest outcomes from treatment  
treat_lower <- treat_sorted %>% slice(1:(n() - n_trim))

# Step 3: Estimate bounds
upper_bound <- mean(treat_upper$earnings_endline, na.rm = TRUE) -
  mean(control_attrited$earnings_endline, na.rm = TRUE)

lower_bound <- mean(treat_lower$earnings_endline, na.rm = TRUE) -
  mean(control_attrited$earnings_endline, na.rm = TRUE)

cat("\n=== Lee Bounds ===\n")
cat("Lower bound:", round(lower_bound, 0), "Ksh\n")
cat("Upper bound:", round(upper_bound, 0), "Ksh\n")
cat("Attrited ITT:", round(attrited_coef, 0), "Ksh\n")
cat("Original ITT (no attrition):", round(coef(itt_simple)["offered"], 0), "Ksh\n")


# Q12: Cost-Benefit Analysis

# Parameters given:
# Cost per participant: $1,500
# Exchange rate: $1 = 79 Ksh
# Time horizon: 2 years (24 months)
# Monthly discount rate: 1%
# Use the LATE (ToT) as the benefit estimate — effect on participants

# We use the IV ToT estimate (4,591 Ksh/month) as our best estimate
# of the causal effect on compliers (actual participants)

cost_per_participant_usd <- 1500
exchange_rate            <- 79       # Ksh per dollar
cost_per_participant_ksh <- cost_per_participant_usd * exchange_rate

monthly_discount_rate    <- 0.01
n_months                 <- 24
tot_estimate_ksh         <- coef(iv_simple)["treated"]  # 4,591 Ksh/month

cat("=== Cost-Benefit Analysis ===\n\n")
cat("Cost per participant: $1,500 =", cost_per_participant_ksh, "Ksh\n")
cat("Monthly benefit (ToT estimate):", round(tot_estimate_ksh, 0), "Ksh\n")
cat("Discount rate:", monthly_discount_rate * 100, "% per month\n")
cat("Time horizon:", n_months, "months\n\n")

# --- Present Value of Benefits ---
# PV = sum of monthly_benefit / (1 + r)^t for t = 1 to 24
months      <- 1:n_months
pv_benefits <- sum(tot_estimate_ksh / (1 + monthly_discount_rate)^months)

cat("PV of benefits over 24 months:", round(pv_benefits, 0), "Ksh\n")
cat("PV of benefits in USD:", round(pv_benefits / exchange_rate, 0), "$\n\n")

# --- Net Present Value ---
npv <- pv_benefits - cost_per_participant_ksh
cat("Cost:", cost_per_participant_ksh, "Ksh\n")
cat("NPV:", round(npv, 0), "Ksh\n")
cat("NPV in USD:", round(npv / exchange_rate, 0), "$\n\n")

# --- Benefit-Cost Ratio ---
bcr <- pv_benefits / cost_per_participant_ksh
cat("Benefit-Cost Ratio:", round(bcr, 3), "\n\n")

# --- Sensitivity check: use lower and upper bounds ---
# Lower: use ITT instead of ToT (conservative)
# Upper: use matching ToT (7,588 Ksh)
itt_ksh     <- coef(itt_simple)["offered"]
match_ksh   <- coef(tot_match)["treated"]

pv_itt     <- sum(itt_ksh   / (1 + monthly_discount_rate)^months)
pv_match   <- sum(match_ksh / (1 + monthly_discount_rate)^months)

npv_itt    <- pv_itt   - cost_per_participant_ksh
npv_match  <- pv_match - cost_per_participant_ksh

bcr_itt    <- pv_itt   / cost_per_participant_ksh
bcr_match  <- pv_match / cost_per_participant_ksh

cat("=== Sensitivity Analysis ===\n")
cat(sprintf("%-25s %10s %10s %8s\n", "Scenario", "PV Ben (Ksh)", "NPV (Ksh)", "BCR"))
cat(sprintf("%-25s %10.0f %10.0f %8.3f\n", 
            "ITT (conservative)",    pv_itt,   npv_itt,   bcr_itt))
cat(sprintf("%-25s %10.0f %10.0f %8.3f\n", 
            "IV ToT (preferred)",    pv_benefits, npv,    bcr))
cat(sprintf("%-25s %10.0f %10.0f %8.3f\n", 
            "Matching ToT (upper)",  pv_match, npv_match, bcr_match))
