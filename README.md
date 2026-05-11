# rct-analysis-kenya-enterprise
R Code for one of my assignments from Designing Field Experiments.
# Cluster-Randomized Trial Analysis: Kenya Enterprise Program
## GPEC 464 — Designing Field Experiments | UC San Diego GPS

## Overview
This project analyzes a two-round cluster-randomized trial of an entrepreneurship 
training program in rural Kenya. The analysis covers the full lifecycle of an RCT 
evaluation: power analysis at the design stage, balance testing, ITT and ToT estimation 
via IV and matching, randomization inference, heterogeneous treatment effects, Lee 
Bounds under differential attrition, and cost-benefit analysis.

## Research Question
Does a village-level entrepreneurship training program increase monthly earnings for 
participants in rural Kenya?

## Data
- **Design:** Cluster-randomized trial, treatment assigned at the village level
- **Sample:** 804 individuals across 168 villages, 2 rounds (baseline + endline)
- **Outcome:** Monthly earnings in Kenyan Shillings (Ksh)
- **Treatment:** Binary — offered training (ITT) / received training (ToT)
- **Covariates:** Baseline earnings, trusting, risk amount, age, urban status, 
  highest grade completed
- **Compliance:** One-sided non-compliance; uptake rate = 38.6% among offered

## Methods & Key Results

### Power Analysis (Pre-Treatment)
- J = 168 villages, n ≈ 4.79 per cluster, ICC = 0.056
- MDE at 80% power: 1,028 Ksh/month without covariates; 982 Ksh/month with covariates
- Covariate adjustment (R² = 0.088) provides only modest precision gain

### Balance Testing
- Two-sample t-tests on six baseline covariates
- No covariate significant at 5% level — randomization is balanced

### ITT Estimation
- Simple OLS: 1,774 Ksh/month (p = 0.133)
- ANCOVA (with baseline covariates): 1,744 Ksh/month (p = 0.121)
- Insignificance driven by low uptake diluting the ITT across compliers and never-takers

### Treatment on the Treated (ToT)
- **IV (2SLS):** 4,591 Ksh/month — instruments "treated" with "offered"; 
  recovers LATE for compliers; not statistically significant
- **Matching:** 7,588 Ksh/month (p = 0.012) — 1-to-1 nearest neighbor PSM 
  with replacement; likely upward-biased due to selection on unobservables
- Gap between IV and matching estimates (>50%) indicates positive selection 
  among compliers

### Randomization Inference
- 500 village-level permutations
- RI SE (1,351 Ksh) converges with clustered parametric SE (1,357 Ksh), 
  validating the inference approach

### Heterogeneous Treatment Effects
- Six subgroup interactions tested; only urban approaches significance (p = 0.064)
- After Bonferroni correction: no heterogeneous effects survive multiple testing
- Conclusion: ITT appears relatively uniform across subgroups

### Lee Bounds (Differential Attrition)
- Imposed 20% treatment / 40% control attrition artificially
- Trimming fraction: 25%; bounds range from -2,795 to +2,231 Ksh
- Bounds include zero — too wide to be policy-informative under differential attrition

### Cost-Benefit Analysis
| Scenario | PV Benefits | NPV | BCR |
|---|---|---|---|
| ITT (lower bound) | 37,685 Ksh | -80,815 Ksh | 0.32 |
| IV ToT (preferred) | 97,538 Ksh | -20,962 Ksh | 0.82 |
| Matching ToT (upper) | 161,190 Ksh | +42,690 Ksh | 1.36 |

Program does not pass a cost-benefit test under the preferred IV estimate (BCR = 0.82). 
Conclusion is sensitive to estimator choice; scaling not recommended without stronger evidence.

## Files
- `dfe_hw1.R` — Full analysis script covering Q1–Q13
- `Q3b_balance_table.html` — Balance of randomization on baseline covariates
- `Q5_ITT_table.html` — ITT estimates: Simple OLS vs. ANCOVA
- `Q7_ToT_IV_table.html` — ToT estimates: IV (2SLS)
- `Q8_matching_table.html` — ToT estimates: Matching estimator
- `Q10_heterogeneous_ITT.html` — Heterogeneous ITT with multiple inference correction
- `power_curve_HW1.png` — Power curve under two covariate scenarios
- `RI_distribution_HW1.png` — Randomization inference null distribution

## R Packages
```r
tidyverse, WebPower, ICC, fishmethods, stargazer, modelsummary,
AER, sandwich, lmtest, margins, kableExtra, MatchIt
```

## How to Run
1. Place `HW1.csv` in your working directory
2. Update `setwd()` path in the script header
3. Run sequentially — all output tables and figures save to the working directory
