library(dplyr)
library(tableone)
library(xtable)
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
  mutate(outcome = readmcount30) %>%
  rename("treatment" = "intervention") %>% 
  rename("period" = "time") %>% 
  rename("cluster" = "id_site") %>% 
  rename("D" = "noted_consult") 

# Extract relevant covariates and the treatment.
# intervention is the treatment
# Extract smaller set of relevant covariates and the treatment based on PI input
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

# bind the one_hot with the outcome, treatment
data <- data.frame(D = clean_data[,"D"], outcome = clean_data[,"outcome"] ,
                   treatment = clean_data[,"treatment"], period = clean_data[,"period"], cluster = clean_data[,"cluster"])
data$treatment <- data$treatment - 1
data <- bind_cols(data, demean_one_hot_covariate_matrix)

covariate_names = colnames(data)[6:ncol(data)]

# ANCOVAI and ANCOVA III and unadjusted
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
ANCOVAI.point.est.confirm <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.0234)
ANCOVAI.point.est.confirm$CR3_deviate
ANCOVAIII.point.est.confirm <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.02)
ANCOVAIII.point.est.confirm$CR3_deviate
unANCOVA.point.est.confirm <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = 0.0078125)
unANCOVA.point.est.confirm$CR3_deviate

# p-value for test at 0, one-sided
ANCOVAI.null <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = 0)
ANCOVAIII.null <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = 0)
unANCOVA.null <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = 0)

pt(ANCOVAI.ITT$CR3_deviate, df = 11 - 2)
pt(ANCOVAIII.ITT$CR3_deviate, df = 11 - 2)
pt(unANCOVA.ITT$CR3_deviate, df = 11 - 2)

# ANCOVA CIs
threshold <- qt(1 - alpha / 2, df = 11 - 2)
threshold_90 <- qt(1 - 0.1 / 2, df = 11 - 2)

# ANCOVAI
# 90%
ANCOVAI.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = 0.0285)
ANCOVAI.upper.ci.CR3.90$CR3_deviate
ANCOVAI.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.109)
ANCOVAI.lower.ci.CR3.90$CR3_deviate
# 95%
ANCOVAI.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = 0.0427)
ANCOVAI.upper.ci.CR3.95$CR3_deviate
ANCOVAI.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAI, lambda = -0.157)
ANCOVAI.lower.ci.CR3.95$CR3_deviate

# ANCOVAIII
# 90%
ANCOVAIII.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = 0.026)
ANCOVAIII.upper.ci.CR3.90$CR3_deviate
ANCOVAIII.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.093)
ANCOVAIII.lower.ci.CR3.90$CR3_deviate

# 95%
ANCOVAIII.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = 0.0394)
ANCOVAIII.upper.ci.CR3.95$CR3_deviate
ANCOVAIII.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = ANCOVAIII, lambda = -0.132)
ANCOVAIII.lower.ci.CR3.95$CR3_deviate

# unadjusted ANCOVA
# 90%
unANCOVA.upper.ci.CR3.90 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.0088)
unANCOVA.upper.ci.CR3.90$CR3_deviate
unANCOVA.lower.ci.CR3.90 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.16)
unANCOVA.lower.ci.CR3.90$CR3_deviate

# 95%
unANCOVA.upper.ci.CR3.95 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = 0.007)
unANCOVA.upper.ci.CR3.95$CR3_deviate
unANCOVA.lower.ci.CR3.95 <- effect_ratio_test(data = data, test_function = unadjusted_ANCOVA, lambda = -0.22)
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


