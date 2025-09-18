library(dplyr)
library(reshape2)
library(tidyr)
library(xtable)

alpha = 0.05
deviate_threshold = qnorm(1 - alpha / 2)
qt(1 - alpha / 2, df = 10)
folds = 10

cluster_periods <- c("I11J10", "I12J5", "I30J5", "I60J5", "I90J5")
HT_methods <- c("unadj_HT", "reg_adj_HT", "lmm_adj_HT", "reg_adj_full_HT") 
ANCOVA_methods <- c("unadj_ANCOVA", "ANCOVAI", "ANCOVAIII")
ivmodel_methods <- c("unadj_ivmodel", "adj_ivmodel")
inf_vector <- c("TRUE", "FALSE")

column_names <- c("cluster-period", "metric", "iteration", "value")

# load in HT results
HT_list <- list()
for (cluster_period in cluster_periods) {
  for (method in HT_methods) {
    for (informative in inf_vector) {
      sim_name <- paste(cluster_period, method, informative, sep = "_")
      results_dir <- paste0("simulations/results/", sim_name)
      current_rds <- readRDS(paste0(results_dir, "/", sim_name, "_results.rds"))
      current_df <- melt(current_rds)
      colnames(current_df) <- column_names
      current_df$informative <- informative
      current_df$method <- method
      HT_list[[length(HT_list) + 1]] <- current_df
      
    }
  }
}
HT_df <- do.call(rbind, HT_list) |> pivot_wider(names_from = metric, values_from = value)


# process function
summarise_results_HT <- function(data) {
  data <- data |> mutate(typeI = abs(deviate) > deviate_threshold,
                         cons.typeI = abs(cons.deviate) > deviate_threshold,
                         power = abs(zero.deviate) > deviate_threshold,
                         cons.power = abs(cons.zero.deviate) > deviate_threshold,
                         bias = point.est - true.lambda,
                         squared_bias = (point.est - true.lambda)^2 )
  summary_data <- data |> group_by(`cluster-period`, method, informative) |> summarize(mean.bias = mean(bias),
                                                                  mse = mean(squared_bias),
                                                                  type.I.rate = ifelse(sum(!is.na(typeI)) > 0, mean(typeI, na.rm = TRUE), NA),
                                                                  proportion.na.typeI = sum(is.na(typeI)) / n(), 
                                                                  cons.type.I.rate = mean(cons.typeI),
                                                                  power.rate = ifelse(sum(!is.na(power)) > 0, mean(power, na.rm = TRUE), NA),
                                                                  proportion.na.power = sum(is.na(power)) / n(), 
                                                                  cons.power.rate = mean(cons.power))
  return(list(data = data, summary_data = summary_data))
}

HT_object <- summarise_results_HT(HT_df)
HT_results <- HT_object$data
HT_summary <- HT_object$summary_data

HT_summary <- HT_summary %>%
  filter(method != "lmm_adj_HT") %>%
  arrange(desc(informative))

# make latex table
print(xtable(HT_summary, digits = c(0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3)), type = "latex",
      include.rownames = FALSE)

# load in ANCOVA results
ANCOVA_list <- list()
for (cluster_period in cluster_periods) {
  for (method in ANCOVA_methods) {
    for (informative in inf_vector) {
      sim_name <- paste(cluster_period, method, informative, sep = "_")
      results_dir <- paste0("simulations/results/", sim_name)
      current_rds <- readRDS(paste0(results_dir, "/", sim_name, "_results.rds"))
      current_df <- melt(current_rds)
      colnames(current_df) <- column_names
      current_df$informative <- informative
      current_df$method <- method
      current_df$clusters <- as.numeric(sub(".*I(\\d+)J.*", "\\1", cluster_period))
      ANCOVA_list[[length(ANCOVA_list) + 1]] <- current_df
      
    }
  }
}
ANCOVA_df <- do.call(rbind, ANCOVA_list) |> pivot_wider(names_from = metric, values_from = value)

# process function
summarise_results_ANCOVA <- function(data) {
  data <- data |> mutate(CR0.typeI = abs(CR0.deviate) > deviate_threshold,
                         CR3.typeI = abs(CR3.deviate) > deviate_threshold,
                         DB.typeI = abs(DB.deviate) > deviate_threshold,
                         CR3.typeI.t = abs(CR3.deviate) > qt(1 - alpha / 2, df = clusters - 2),
                         CR0.power = abs(CR0.zero.deviate) > deviate_threshold,
                         CR3.power = abs(CR3.zero.deviate) > deviate_threshold,
                         CR3.power.t = abs(CR3.zero.deviate) > qt(1 - alpha / 2, df = clusters - 2),
                         DB.power = abs(DB.zero.deviate) > deviate_threshold,
                         bias = point.est - true.lambda,
                         squared_bias = (point.est - true.lambda)^2 )
  summary_data <- data |> group_by(`cluster-period`, method, informative) |> summarize(mean.bias = mean(bias),
                                                                                       mse = mean(squared_bias),
                                                                                       CR0.type.I.rate = mean(CR0.typeI),
                                                                                       CR3.type.I.rate = ifelse(sum(!is.na(CR3.typeI)) > 0, mean(CR3.typeI, na.rm = TRUE), NA),
                                                                                       proportion.na.CR3.typeI = sum(is.na(CR3.typeI)) / n(), 
                                                                                       DB.type.I.rate = mean(DB.typeI),
                                                                                       CR0.power.rate = mean(CR0.power),
                                                                                       CR3.power.rate = ifelse(sum(!is.na(CR3.power)) > 0, mean(CR3.power, na.rm = TRUE), NA),
                                                                                       proportion.na.CR3.power = sum(is.na(CR3.power)) / n(), 
                                                                                       DB.power.rate = mean(DB.power))
  return(list(data = data, summary_data = summary_data))
}

# process function
summarise_results_ANCOVA_no_na <- function(data) {
  data <- data |> mutate(CR0.typeI = abs(CR0.deviate) > deviate_threshold,
                         CR3.typeI = abs(CR3.deviate) > deviate_threshold,
                         DB.typeI = abs(DB.deviate) > deviate_threshold,
                         CR3.typeI.t = ifelse(is.na(abs(CR3.deviate) > qt(1 - alpha / 2, df = clusters - 2)), 0, 
                                              abs(CR3.deviate) > qt(1 - alpha / 2, df = clusters - 2)),
                         CR0.power = abs(CR0.zero.deviate) > deviate_threshold,
                         CR3.power = abs(CR3.zero.deviate) > deviate_threshold,
                         CR3.power.t = ifelse(is.na(abs(CR3.zero.deviate) > qt(1 - alpha / 2, df = clusters - 2)), 0, 
                                              abs(CR3.zero.deviate) > qt(1 - alpha / 2, df = clusters - 2)),
                         DB.power = abs(DB.zero.deviate) > deviate_threshold,
                         bias = point.est - true.lambda,
                         squared_bias = (point.est - true.lambda)^2 )
  summary_data <- data |> group_by(`cluster-period`, method, informative) |> summarize(mean.bias = mean(bias),
                                                                                       mse = mean(squared_bias),
                                                                                       CR0.type.I.rate = mean(CR0.typeI),
                                                                                       CR3.type.I.rate = mean(CR3.typeI, na.rm = TRUE),
                                                                                       CR3.type.I.t.rate = mean(CR3.typeI.t, na.rm = TRUE),
                                                                                       DB.type.I.rate = mean(DB.typeI),
                                                                                       CR0.power.rate = mean(CR0.power),
                                                                                       CR3.power.rate = mean(CR3.power, na.rm = TRUE),
                                                                                       CR3.power.t.rate = mean(CR3.power.t, na.rm = TRUE),
                                                                                       DB.power.rate = mean(DB.power))
  return(list(data = data, summary_data = summary_data))
}
ANCOVA_object <- summarise_results_ANCOVA_no_na(ANCOVA_df)
ANCOVA_object <- summarise_results_ANCOVA(ANCOVA_df)
ANCOVA_results <- ANCOVA_object$data
ANCOVA_summary <- ANCOVA_object$summary_data

ANCOVA_summary <- ANCOVA_summary %>%
  arrange(desc(informative))


# make latex table
print(xtable(ANCOVA_summary, digits = c(0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)), type = "latex",
      include.rownames = FALSE) 

# load in HT results
ivmodel_list <- list()
for (cluster_period in cluster_periods) {
  for (method in ivmodel_methods) {
    for (informative in inf_vector) {
      dummy_df <- data.frame()
      for (j in 1:folds) {
        sim_name <- paste(cluster_period, method, informative, j, sep = "_")
        results_dir <- paste0("simulations/results/", sim_name)
        current_rds <- readRDS(paste0(results_dir, "/", sim_name, "_results.rds"))
        current_df <- melt(current_rds)
        current_df$Dimension3 <- current_df$Dimension3 + (j - 1) * 100
        colnames(current_df) <- column_names
        current_df$informative <- informative
        current_df$method <- method
        dummy_df <- rbind(dummy_df,current_df)
  
      }
      ivmodel_list[[length(ivmodel_list) + 1]] <- dummy_df
    }
  }
}
ivmodel_df <- do.call(rbind, ivmodel_list) |> pivot_wider(names_from = metric, values_from = value)


# process function
summarise_results_ivmodel <- function(data) {
  data <- data |> mutate(CLR.typeI = (true.lambda < lower.ci.CLR) | (true.lambda > upper.ci.CLR),
                         CLR.power = (0 < lower.ci.CLR) | (0 > upper.ci.CLR),
                         AR.typeI = (true.lambda < lower.ci.AR) | (true.lambda > upper.ci.AR),
                         AR.power = (0 < lower.ci.AR) | (0 > upper.ci.AR),
                         LIML.typeI = (true.lambda < lower.ci.LIML) | (true.lambda > upper.ci.LIML),
                         LIML.power = (0 < lower.ci.LIML) | (0 > upper.ci.LIML),
                         LIML.bias = point.est.LIML - true.lambda,
                         LIML.squared_bias = (point.est.LIML - true.lambda)^2,
                         Fuller.typeI = (true.lambda < lower.ci.Fuller) | (true.lambda > upper.ci.Fuller),
                         Fuller.power = (0 < lower.ci.Fuller) | (0 > upper.ci.Fuller),
                         Fuller.bias = point.est.Fuller - true.lambda,
                         Fuller.squared_bias = (point.est.Fuller - true.lambda)^2)
  summary_data <- data |> group_by(`cluster-period`, method, informative) |> summarize(LIML.mean.bias = mean(LIML.bias),
                                                                                       LIML.mse = mean(LIML.squared_bias),
                                                                                       Fuller.mean.bias = mean(Fuller.bias),
                                                                                       Fuller.mse = mean(Fuller.squared_bias),
                                                                                       CLR.type.I.rate = mean(CLR.typeI),
                                                                                       AR.type.I.rate = mean(AR.typeI),
                                                                                       LIML.type.I.rate = mean(LIML.typeI),
                                                                                       Fuller.type.I.rate = mean(Fuller.typeI),
                                                                                       CLR.power.rate = mean(CLR.power),
                                                                                       AR.power.rate = mean(AR.power),
                                                                                       LIML.power.rate = mean(LIML.power),
                                                                                       Fuller.power.rate = mean(Fuller.power))
  return(list(data = data, summary_data = summary_data))
}

ivmodel_object <- summarise_results_ivmodel(ivmodel_df)
ivmodel_results <- ivmodel_object$data
ivmodel_summary <- ivmodel_object$summary_data


ivmodel_summary <- ivmodel_summary %>%
  arrange(desc(informative)) %>% select(-c(LIML.mean.bias, LIML.mse, AR.type.I.rate,
                                           LIML.type.I.rate, AR.power.rate,
                                           LIML.power.rate))

# make latex table
print(xtable(ivmodel_summary, digits = c(0, 0, 0, 0, 3, 3, 3, 3, 3, 3)), type = "latex",
      include.rownames = FALSE) 
