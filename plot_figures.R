library(dplyr)
library(reshape2)
library(tidyr)
library(xtable)
library(ggplot2)

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
                                                                                       proportion.na.typeI = sum(is.nan(deviate)) / n(), 
                                                                                       cons.type.I.rate = mean(cons.typeI),
                                                                                       power.rate = ifelse(sum(!is.na(power)) > 0, mean(power, na.rm = TRUE), NA),
                                                                                       proportion.na.power = sum(is.nan(zero.deviate)) / n(), 
                                                                                       cons.power.rate = mean(cons.power))
  return(list(data = data, summary_data = summary_data))
}

HT_object <- summarise_results_HT(HT_df)
HT_results <- HT_object$data
HT_summary <- HT_object$summary_data

HT_summary <- HT_summary %>%
  filter(method != "lmm_adj_HT") %>%
  filter(informative == TRUE)

HT_summary <- HT_summary |> pivot_longer(
  cols = ends_with("type.I.rate") | ends_with("power.rate"),
  names_to = c("scenario", ".value"),
  names_pattern = "(cons)?\\.?(.+)"
) |> mutate(scenario = ifelse(scenario == "cons", "cons", "base")) |>
                           mutate(clusters = as.numeric(substring(`cluster-period`, 2, 3))) |> 
                           rename(type.I = type.I.rate, power = power.rate) |>
  pivot_longer(cols = c(mse, type.I, power), names_to = "metric") |>
                           mutate(cons = scenario == "cons") |>
                           mutate(method =  case_when(
                                  method %in% c("reg_adj_HT", "reg_adj_HT_cons")  ~ "reg_adj",
                                  method %in% c("reg_adj_full_HT", "reg_adj_full_HT_cons")  ~ "reg_adj_full",
                                  method %in% c("unadj_HT", "unadj_HT_cons")  ~ "unadj",
                                  TRUE ~ method))
HT_plot <- ggplot(HT_summary |> drop_na() |> mutate(method_cons = paste0(method,"_",cons)) |>
                    filter(proportion.na.typeI < 0.99 | cons == TRUE | metric == "mse"),
                  aes(x = `cluster-period`, y = value, color = method, linetype = cons,
                      group = method_cons)) +
  geom_line() + 
  ggtitle("Performance of Horvitz-Thompson methods across no. clusters with informative cluster size") +
  facet_wrap(~metric, ncol = 3, scales = "free_y") + theme_bw() 
HT_plot
ggsave("figures/ht-informative.pdf", plot = HT_plot, width = 8, height = 3)

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
                                                                                       CR3.t.type.I.rate = mean(CR3.typeI.t, na.rm = TRUE),
                                                                                       DB.type.I.rate = mean(DB.typeI),
                                                                                       CR0.power.rate = mean(CR0.power),
                                                                                       CR3.power.rate = mean(CR3.power, na.rm = TRUE),
                                                                                       CR3.t.power.rate = mean(CR3.power.t, na.rm = TRUE),
                                                                                       DB.power.rate = mean(DB.power))
  return(list(data = data, summary_data = summary_data))
}
ANCOVA_object <- summarise_results_ANCOVA_no_na(ANCOVA_df)
# ANCOVA_object <- summarise_results_ANCOVA(ANCOVA_df)
ANCOVA_results <- ANCOVA_object$data
ANCOVA_summary <- ANCOVA_object$summary_data

ANCOVA_summary <- ANCOVA_summary %>%
  filter(informative == TRUE)


ANCOVA_summary <- ANCOVA_summary |> pivot_longer(
  cols = ends_with(".type.I.rate") | ends_with(".power.rate"),
  names_to = c("meth", ".value"),
  names_pattern = "(.*)\\.(type\\.I\\.rate|power\\.rate)"
) |> mutate(clusters = as.numeric(substring(`cluster-period`, 2, 3))) |> 
  rename(type.I = type.I.rate,
         power = power.rate,
         var = meth) |>
  pivot_longer(cols = c(mse, type.I, power), names_to = "metric") |>
  mutate(method =  case_when(
    method == "ANCOVAI"  ~ "ANCI",
    method == "ANCOVAIII"  ~ "ANCIII",
    method == "unadj_ANCOVA"  ~ "unadj",
    TRUE ~ method))
ANCOVA_plot <- ggplot(ANCOVA_summary |> drop_na() |> filter(var != "CR3") |>
                        mutate(method_var = paste0(method,"_",var)),
                  aes(x = `cluster-period`, y = value, color = method, linetype = var,
                      group = method_var)) +
  geom_line() + 
  ggtitle("Performance of ANCOVA methods across no. clusters with informative cluster size") +
  facet_wrap(~metric, ncol = 3, scales = "free_y") + theme_bw() 
ANCOVA_plot
ggsave("figures/ancova-informative.pdf", plot = ANCOVA_plot, width = 8, height = 3)

# load in IV results
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
                                           LIML.power.rate)) %>%
  filter(informative == TRUE)

ivmodel_summary <- ivmodel_summary |> pivot_longer(
  cols = ends_with(".type.I.rate") | ends_with(".power.rate"),
  names_to = c("meth", ".value"),
  names_pattern = "(.*)\\.(type\\.I\\.rate|power\\.rate)"
) |> mutate(clusters = as.numeric(substring(`cluster-period`, 2, 3))) |> 
  rename(type.I = type.I.rate,
         power = power.rate,
         var = meth,
         mse = Fuller.mse) |>
  pivot_longer(cols = c(mse, type.I, power), names_to = "metric") |>
  mutate(method =  case_when(
    method == "adj_ivmodel"  ~ "adj",
    method == "unadj_ivmodel"  ~ "unadj",
    TRUE ~ method))
ivmodel_plot <- ggplot(ivmodel_summary |> drop_na() |> mutate(method_var = paste0(method,"_",var)),
                       aes(x = `cluster-period`, y = value, color = method, linetype = var,
                           group = method_var)) +
  geom_line() + 
  ggtitle("Performance of ivmodel methods across no. clusters with informative cluster size") +
  facet_wrap(~metric, ncol = 3, scales = "free_y") + theme_bw() 
ivmodel_plot
ggsave("figures/ivmodel-informative.pdf", plot = ivmodel_plot, width = 8, height = 3)

# Across method plots
summary <- bind_rows(HT_summary |> mutate(method = paste0(method,"_",cons)), 
                     ANCOVA_summary |> mutate(method = paste0(method,"_",var)), 
                     ivmodel_summary |> mutate(method = paste0(method,"_",var)))
# keep a subset
subset_names <- c("reg_adj_FALSE", "reg_adj_TRUE", "reg_adj_full_TRUE", "ANCI_CR3.t",
                  "ANCIII_CR3.t", "adj_Fuller")
subset_summary <- summary |> filter(method %in% subset_names)
summary_plot <- ggplot(subset_summary,
                       aes(x = `cluster-period`, y = value, color = method, 
                           group = method)) +
  geom_line() + 
  ggtitle("Performance of best methods across no. clusters with informative cluster size") +
  facet_wrap(~metric, ncol = 3, scales = "free_y") + theme_bw() 
summary_plot
ggsave("figures/summary-informative.pdf", plot = summary_plot, width = 8, height = 3)
