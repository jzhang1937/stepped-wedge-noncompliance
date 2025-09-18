library(dplyr)
library(tableone)
library(xtable)
library(ggplot2)
library(purrr)
library(tidyr)
library(broom)
library(sandwich)
library(lmtest)
library(tibble)
source("estimators.R")
source("estimator_helpers.R")

redat= read.csv('AnalysisFile.csv')

### mITT sample analysis
re_mitt=redat[redat$mITT==1,]


### ranking death at 99th percentile ~ 26 days
re_mitt$ranked_los= ifelse(re_mitt$new_death_dispo_add==1,quantile(re_mitt$outcome_los_hrs_add, probs = 0.99),re_mitt$outcome_los_hrs_add)
re_mitt$time= as.factor(re_mitt$step)
re_mitt$intervention=relevel(as.factor(re_mitt$intervention), ref='0')

# process data
data <- re_mitt %>% mutate(time = as.numeric(time) - 1) %>%
  mutate(outcome = log(ranked_los)) %>%
  rename("treatment" = "intervention") %>% 
  rename("period" = "time") %>% 
  rename("cluster" = "id_site") %>% 
  rename("D" = "noted_consult") 

# Extract relevant covariates and the treatment.
# intervention is the treatment
relevant <- c("treatment","D","outcome","cluster","age","mari_cat","admi_cat","elixscore_ahrq","sample_dem","sample_esrd","sample_copd","mITT_TimeSince1st","enrolled_icu", "period")

clean_data <- cbind(data[,relevant])
clean_data <- na.omit(clean_data)
clean_data$treatment <- as.numeric(clean_data$treatment)

covariates <- c("age","admi_cat","elixscore_ahrq","sample_dem","sample_esrd","sample_copd","mITT_TimeSince1st","enrolled_icu")

# Create one hot-encoded covariate matrix.
covariate_matrix <- clean_data[,covariates]
covariate_matrix$admi_cat <- as.factor(covariate_matrix$admi_cat)

# Convert factor variables to one-hot encoding
one_hot_covariate_matrix <- model.matrix(~ . - 1, data = covariate_matrix)
# get rid of redundant categorical column to prevent collinear
one_hot_covariate_matrix <- as.matrix(one_hot_covariate_matrix[,-which(colnames(one_hot_covariate_matrix) == "admi_cat1")])

# Get rid of period mean
df <- data.frame(cbind(one_hot_covariate_matrix, period = factor(clean_data$period)))
demean_one_hot_covariate_matrix <- df %>% mutate(period = factor(period)) %>%
  group_by(period) %>%
  mutate_at(vars(-period), funs(. - mean(.)))
demean_one_hot_covariate_matrix = as.matrix(demean_one_hot_covariate_matrix[,-which(colnames(demean_one_hot_covariate_matrix) == "period")])

# bind the one_hot with the outcome, treatment.
data <- data.frame(D = clean_data[,"D"], outcome = clean_data[,"outcome"] ,
                   treatment = clean_data[,"treatment"], period = clean_data[,"period"], cluster = clean_data[,"cluster"])
data$treatment <- data$treatment - 1
data <- bind_cols(data, demean_one_hot_covariate_matrix)

covariate_names = colnames(data)[6:ncol(data)]

# Make a histogram of cluster-period sizes in the rollout periods
cluster_period_counts <- rename(count(data |> filter(period != 0 & period != 11), period, cluster), `cluster-period size` = n)
time_period_plot <- ggplot(data = cluster_period_counts, aes(x = `cluster-period size`)) +
  geom_histogram(bins = 50)+theme_bw(base_size = 18)
print(time_period_plot)
ggsave("analysis/cluster_period_counts.pdf", plot = time_period_plot)

# Make a table 1 excluding pre and post rollout
dummy_data <- data.frame(D = clean_data[,"D"], outcome = clean_data[,"outcome"] ,
                         treatment = clean_data[,"treatment"], period = clean_data[,"period"], cluster = clean_data[,"cluster"])
dummy_data$treatment <- dummy_data$treatment - 1
data_table1 <- bind_cols(dummy_data, df[,1:10]) %>% filter(period != 0 & period != 11)

table <- CreateTableOne(vars = colnames(data_table1)[6:15], strata = "treatment", data = data_table1, 
                        includeNA = FALSE, addOverall = FALSE, test = FALSE)

# Print table
# Convert TableOne to a printable format
table_df <- print(table, showAllLevels = TRUE, format = "p", quote = FALSE, smd = TRUE,
                  noSpaces = TRUE, catDigits = 2, contDigits = 2)

# Convert to data frame for xtable
table_df <- as.data.frame(table_df)

# Convert to LaTeX using xtable
table_latex <- xtable(table_df)

# Print LaTeX output
print(table_latex, type = "latex", include.rownames = TRUE)

# ANCOVAI and ANCOVA III and unadjusted
source("estimators.R")
source("estimator_helpers.R")

alpha = 0.05

# ITT
data_ITT <- data
data_ITT$D <- NULL

ANCOVAI.ITT <- ANCOVAI(data = data_ITT)
ANCOVAIII.ITT <- ANCOVAIII(data = data_ITT)
unANCOVA.ITT <- unadjusted_ANCOVA(data = data_ITT)

# ANCOVA point estimates
ANCOVAI.point.est <- effect_ratio_point_estimate(data = data, test_function = ANCOVAI, 
                                                 tol = 0.01, max_iter = 20, min = -1, max = 1)
ANCOVAIII.point.est <- effect_ratio_point_estimate(data = data, test_function = ANCOVAIII,
                                                   tol = 0.01, max_iter = 20, min = -1, max = 1)
unANCOVA.point.est <- effect_ratio_point_estimate(data = data, test_function = unadjusted_ANCOVA, 
                                                 tol = 0.01, max_iter = 20, min = -1, max = 1)
ANCOVAI.point.est.confirm <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.18)
ANCOVAI.point.est.confirm$CR3_deviate
ANCOVAIII.point.est.confirm <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.175)
ANCOVAIII.point.est.confirm$CR3_deviate
unANCOVA.point.est.confirm <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.148)
unANCOVA.point.est.confirm$CR3_deviate

pt(ANCOVAI.ITT$CR3_deviate, df = 11 - 2)
pt(ANCOVAIII.ITT$CR3_deviate, df = 11 - 2)
pt(unANCOVA.ITT$CR3_deviate, df = 11 - 2)

# ANCOVA CIs
threshold <- qt(1 - alpha / 2, df = 11 - 2)
threshold_90 <- qt(1 - 0.1 / 2, df = 11 - 2)

# ANCOVAI
# 90%
ANCOVAI.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = 0.009)
ANCOVAI.upper.ci.CR3.90$CR3_deviate
ANCOVAI.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.455)
ANCOVAI.lower.ci.CR3.90$CR3_deviate
# 95%
ANCOVAI.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = 0.068)
ANCOVAI.upper.ci.CR3.95$CR3_deviate
ANCOVAI.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.596)
ANCOVAI.lower.ci.CR3.95$CR3_deviate

# ANCOVAIII
# 90%
ANCOVAIII.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = 0.002)
ANCOVAIII.upper.ci.CR3.90$CR3_deviate
ANCOVAIII.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.462)
ANCOVAIII.lower.ci.CR3.90$CR3_deviate

# 95%
ANCOVAIII.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = 0.052)
ANCOVAIII.upper.ci.CR3.95$CR3_deviate
ANCOVAIII.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.62)
ANCOVAIII.lower.ci.CR3.95$CR3_deviate

# unadjusted ANCOVA
# 90%
unANCOVA.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = 0.035)
unANCOVA.upper.ci.CR3.90$CR3_deviate
unANCOVA.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.37)
unANCOVA.lower.ci.CR3.90$CR3_deviate

# 95%
unANCOVA.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = 0.101)
unANCOVA.upper.ci.CR3.95$CR3_deviate
unANCOVA.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.488)
unANCOVA.lower.ci.CR3.95$CR3_deviate


# one-hot encoded period
one_hot_period_matrix <- model.matrix(~ . - 1, data = data.frame(as.factor(data$period)))
# ivmodel adjusted (adjusted for period and covariates)
adjusted.ivmodel.90 <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment,
                               X = cbind(one_hot_period_matrix[,-1], data[,covariate_names]), clusterID = data$cluster, alpha = 0.1)
adjusted.ivmodel.95 <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment,
                               X = cbind(one_hot_period_matrix[,-1], data[,covariate_names]), clusterID = data$cluster, alpha = 0.05)

# ivmodel unadjusted (adjusted for period)
unadjusted.ivmodel.90 <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment, X = one_hot_period_matrix[,-1],
                                 clusterID = data$cluster, alpha = 0.1)

unadjusted.ivmodel.95 <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment, X = one_hot_period_matrix[,-1],
                                 clusterID = data$cluster, alpha = 0.05)

### HEURISTIC CHECK of TREATMENT DURATION IRRELEVANCE ###
# Use raw covariates
data_raw <- data.frame(D = clean_data[,"D"], outcome = clean_data[,"outcome"] ,
                   treatment = clean_data[,"treatment"], period = clean_data[,"period"], cluster = clean_data[,"cluster"])
data_raw$treatment <- data_raw$treatment - 1
data_raw <- bind_cols(data_raw, one_hot_covariate_matrix)
# treated
treated_rollout <- data_raw |> filter(period != 0 & period != 11 & treatment == 1) |>
  group_by(cluster) %>%
  mutate(adopt = min(period),
         time_on_treatment = period - adopt) %>%
  ungroup()

# List of outcomes to test
outcomes <- c("D", "outcome")
treated_outcome_results <- map_dfr(outcomes, function(.outcome) {
  
  # Define model formula
  formula <- as.formula(paste(.outcome, "~ time_on_treatment + age + admi_cat2 +
                               admi_cat31 + admi_cat32 + elixscore_ahrq + sample_dem +
                               sample_esrd + sample_copd + mITT_TimeSince1st + enrolled_icu"))
  
  # Fit model
  model <- if (.outcome == "D") {
    glm(formula, data = treated_rollout, family = binomial)
  } else {
    lm(formula, data = treated_rollout)
  }
  
  # Cluster-robust SEs
  robust_vcov <- clubSandwich::vcovCR(model, cluster = treated_rollout$cluster, type = "CR0")
  robust_model <- clubSandwich::coef_test(model, vcov = robust_vcov)
  
  robust_model %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(grepl("time_on_treatment", term)) %>%
    mutate(outcome = .outcome)
})
saveRDS(treated_outcome_results, "analysis/treated_test.rds")

# control 
control_rollout <- data_raw |> filter(period != 0 & period != 11 & treatment == 0) |>
  group_by(cluster) %>%
  mutate(time_on_control = period) %>%
  ungroup()

control_outcome_results <- map_dfr(outcomes, function(.outcome) {
  
  # Define model formula
  formula <- as.formula(paste(.outcome, "~ time_on_control + age + admi_cat2 +
                               admi_cat31 + admi_cat32 + elixscore_ahrq + sample_dem +
                               sample_esrd + sample_copd + mITT_TimeSince1st + enrolled_icu"))
  
  # Fit model
  model <- if (.outcome == "D") {
    glm(formula, data = control_rollout, family = binomial)
  } else {
    lm(formula, data = control_rollout)
  }
  
  # Cluster-robust SEs
  robust_vcov <- clubSandwich::vcovCR(model, cluster = control_rollout$cluster, type = "CR0")
  robust_model <- clubSandwich::coef_test(model, vcov = robust_vcov)
  
  robust_model %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(grepl("time_on_control", term)) %>%
    mutate(outcome = .outcome)
})
saveRDS(control_outcome_results, "analysis/control_test.rds")

# print summary table
xtable(rbind(treated_test, control_test), digits = 3)
