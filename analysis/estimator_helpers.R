library(clubSandwich)
library(lme4)
library(MASS) # Ensure the MASS package is loaded
library(dplyr)
# helpers for the estimator functions

# (linear) regression adjustment function using pre and post data.
reg_adjust_HT <- function(data) {
  # find pre and post rollout
  pre_rollout <- min(data$period)
  post_rollout <- max(data$period)
  
  # subset based on pre and post.
  data0 <- data %>% dplyr::filter(period == pre_rollout) %>% dplyr::select(-c(period,treatment,cluster)) 
  data1 <- data %>% dplyr::filter(period == post_rollout) %>% dplyr::select(-c(period,treatment,cluster)) 
  # run regression for treated
  model1 <- lm(outcome ~ ., data = data1 )
  
  g1 <- function(new_data) {
    predict(model1, newdata = new_data)
  }
  
  # run regression for control
  model0 <- lm(outcome ~ ., data = data0 )
  
  g0 <- function(new_data) {
    predict(model0, newdata = new_data)
  }
  
  return(list(g1 = g1, g0 = g0))
  
}

# (linear) regression adjustment function using all data
reg_adjust_HT_full <- function(data) {
  data <- data %>% mutate(period = as.factor(period))
  # subset based on treatment
  data0 <- data %>% dplyr::filter(treatment == 0) %>% dplyr::select(-c(treatment,cluster)) 
  data1 <- data %>% dplyr::filter(treatment == 1) %>% dplyr::select(-c(treatment,cluster)) 
  # run regression for treated
  model1 <- lm(outcome ~ ., data = data1 )
  
  g1 <- function(new_data) {
    predict(model1, newdata = new_data)
  }
  
  # run regression for control
  model0 <- lm(outcome ~ ., data = data0 )
  
  g0 <- function(new_data) {
    predict(model0, newdata = new_data)
  }
  
  return(list(g1 = g1, g0 = g0))
  
}

# (linear) mixed model adjustment function using all data
lmm_HT <- function(data) {
  data <- data %>% mutate(period = as.factor(period)) %>% 
    mutate(cluster = as.factor(cluster))
  
  # Get all the column names except for 'outcome' and 'cluster'
  predictors <- setdiff(names(data), c("outcome", "cluster"))
  
  # Create the formula dynamically, adding random effects for 'cluster'
  formula <- as.formula(paste("outcome ~", paste(predictors, collapse = " + "), "+ (1 | cluster)"))
  
  # Fit the mixed-effects model using lmer
  lmm_model <- lmer(formula, data = data)
  
  or <- function(new_data) {
    predict(lmm_model, newdata = new_data, re.form = NA)
  }
  
  return(or)
  
}

# transform the data to prepare it for doing inference on the effect ratio
# using the estimating equation.
effect_ratio_transform <- function(data, lambda = 0) {
  transformed_data <- data
  transformed_data$outcome <- data$outcome - data$D * lambda
  transformed_data$D <- NULL
  return(transformed_data)
}

# test effect ratio at a specific lambda
effect_ratio_test <- function(data, alpha = 0.05, test_function, lambda = 0) {
  transformed_data <- effect_ratio_transform(data = data, lambda = lambda)
  return(test_function(data = transformed_data))
}

effect_ratio_point_estimate <- function(data, min = NA, max = NA,
                                        test_function, tol = 0.001, 
                                        max_iter = 100) {
  if(is.na(min) | is.na(max)) {
    range <- max(data$outcome) - min(data$outcome) / 2
    min <- -range
    max <- range
  }
  
  # Iterate until the desired tolerance is reached or maximum iterations are exceeded
  iter = 0
  while(iter <= max_iter) {
    # Update iteration counter
    iter <- iter + 1
    
    # Evaluate the function at the midpoint
    mid <- (min + max) / 2
    
    if((max-min)/2 < tol) {
      return(mid)
    }
    
    # compute estimate at min
    min_data <- effect_ratio_transform(data = data, lambda = min)
    est_min <- test_function(data = min_data)$est
    
    # compute estimate at max
    max_data <- effect_ratio_transform(data = data, lambda = max)
    est_max <- test_function(data = max_data)$est
    
    # compute estimate at midpoint
    mid_data <- effect_ratio_transform(data = data, lambda = mid)
    est_mid <- test_function(data = mid_data)$est
    
    # Check if the root is on the left or right side
    if (est_min * est_mid < 0) {
      max <- mid
    } else {
      min <- mid
    }
    
    if(iter %% 10 == 0) {
      print(paste("Iteration:", iter, "Midpoint:", mid))
    }
  }
  return(mid)
}

# grid search for effect ratio, two-sided
effect_ratio_grid_search <- function(data, min = -5, max = 5, alpha = 0.05,
                                     test_function, epsilon = 0.01, max_iter = 100) {
  lower <- qnorm(alpha / 2)
  upper <- qnorm(1 - alpha / 2)
  lower_ci <- NA
  upper_ci <- NA
  # grid search for lower
    for (i in 1:max_iter) {
      mid <- (min + max) / 2  # Midpoint
      # transform data at midpoint
      mid_data <- effect_ratio_transform(data = data, lambda = mid)
      # deviate at the midpoint
      deviate_mid <- test_function(data = mid_data)$deviate
      if (abs(deviate_mid - lower) < epsilon || (max - min) / 2 < epsilon) {
        lower_ci <- mid
        break
      }
      # transform data at midpoint
      min_data <- effect_ratio_transform(data = data, lambda = min)
      # deviate at the midpoint
      deviate_min <- test_function(data = min_data)$deviate
      if ((deviate_mid - lower) * (deviate_min - lower) < 0) {
        max <- mid
      } else {
        min <- mid
      }
      if (i %% 10 == 0) {
        paste0("Lower iteration: ", i, ", Midpoint: ", mid)
      }
    }

  # grid search for upper
  for (i in 1:max_iter) {
    mid <- (min + max) / 2  # Midpoint
    # transform data at midpoint
    mid_data <- effect_ratio_transform(data = data, lambda = mid)
    # deviate at the midpoint
    deviate_mid <- test_function(mid_data)$deviate
    if (abs(deviate_mid - upper) < epsilon || (max - min) / 2 < epsilon) {
      upper_ci <- mid
      break
    }
    # transform data at midpoint
    min_data <- effect_ratio_transform(data = data, lambda = min)
    # deviate at the midpoint
    deviate_min <- test_function(min_data)$deviate
    if ((deviate_mid - upper) * (deviate_min - upper) < 0) {
      max <- mid
    } else {
      min <- mid
    }
    if (i %% 10 == 0) {
      print(paste0("Upper iteration: ", i, ", Midpoint: ", mid))
    }
  }
  
  return(c(lower_ci, upper_ci))
  
}

# DB variance for ANCOVA
plug_in_cov <- function(data, residual, weights, I, J, I_tr_seq, I_ct_seq, I_a) {
  
  data$res_x <- residual
  data$weights <- weights
  
  # calculate \bar{u}_j(1) and \bar{u}_j(0)
  u1_list <- lapply(1:J, function(x) weighted.mean(data$res_x[(data$treatment==1) & (data$period==x)],
                                                       data$weights[(data$treatment==1) & (data$period==x)]))
  u0_list <- lapply(1:J, function(x) weighted.mean(data$res_x[(data$treatment==0) & (data$period==x)],
                                                       data$weights[(data$treatment==0) & (data$period==x)]))
  u1 = unlist(u1_list)
  u0 = unlist(u0_list)
  
  wj1 = tapply(data$weights[data$treatment==1], data$period[data$treatment==1], sum) / I_tr_seq
  wj0 = tapply(data$weights[data$treatment==0], data$period[data$treatment==0], sum) / I_ct_seq
  
  # debugging
  # print(u1)
  # print(u0)
  # print(wj1)
  # print(wj0)
  
  cluster_cov_ar <- array(NA, dim = c(J+1, J, J))
  for (a in 1:(J+1)) {
    data_a <- data[data$schedule==a,]
    suppressMessages(data_a_ij <- data_a %>% 
                       group_by(period, cluster) %>% 
                       summarise(mean_res_x = weighted.mean(res_x, weights),
                                 w_ij = sum(weights)))
    uj_mean_a = u1*(a<=c(1:J)) + u0*(a>c(1:J))
    
    # uj_mean_a_list <- lapply(2:(J-1), function(x) weighted.mean(data_a$res_x[data_a$period==x],
    #                                                             data_a$weights[data_a$period==x]))
    # uj_mean_a = unlist(uj_mean_a_list)
    
    # wja = tapply(data_a$weights, data_a$period, sum) / 10
    
    cluster_a_mat_list <- lapply(unique(data_a_ij$cluster), function(x) {
      data_a_ij_tmp <- data_a_ij[data_a_ij$cluster==x,]
      vec <- data_a_ij_tmp$mean_res_x - uj_mean_a
      weight_mat <- diag(((data_a_ij_tmp$w_ij*I_a)/(wj1*I_tr_seq))*(a<=c(1:J)) -
                           ((data_a_ij_tmp$w_ij*I_a)/(wj0*I_ct_seq))*(a>c(1:J)))
      c(vec %*% weight_mat) %o% c(vec %*% weight_mat)
    })
    # If one cluster at a time, 
    if (I_a - 1 == 0) {
      cluster_a_cov <- (1/I_a)*1/(I_a)*Reduce('+', cluster_a_mat_list)
    } else {
      cluster_a_cov <- (1/I_a)*1/(I_a-1)*Reduce('+', cluster_a_mat_list)
    }
    cluster_cov_ar[a,,] <- cluster_a_cov

  }
  cov_mat <- apply(cluster_cov_ar, c(2,3), sum)
  # varpi <- tapply(data$weights, data$period, sum) / sum(data$weights)
  # var_out <- c(varpi %*% cov_mat %*% varpi)
  
  return(cov_mat)
}


# psuedo inverse version of covariance estimation
vcovCR.pseudo.lm <- function (obj, cluster, type, target = NULL, inverse_var = NULL, 
          form = "sandwich", ...) 
{
  if (missing(cluster)) 
    stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) 
    inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR_pseudo(obj, cluster = cluster, type = type, target = target, 
          inverse_var = inverse_var, form = form)
}

vcov_CR_pseudo <- function (obj, cluster, type, target = NULL, inverse_var = FALSE, 
          form = "sandwich", ignore_FE = FALSE) 
{
  cluster <- droplevels(as.factor(cluster))
  alias <- is.na(clubSandwich:::coef_CS(obj))
  X <- clubSandwich:::model_matrix(obj)
  if (any(alias)) {
    X <- X[, !alias, drop = FALSE]
  }
  p <- NCOL(X)
  N <- NROW(X)
  cluster_length <- length(cluster)

  if (cluster_length != N) {
    cluster <- droplevels(handle_vectors(cluster, obj))
    if (length(cluster) != N) {
      stop("Clustering variable must have length equal to the number of rows in the data used to fit obj.")
    }
  }
  if (any(is.na(cluster))) 
    stop("Clustering variable cannot have missing values.")
  J <- nlevels(cluster)
  if (J < 2) 
    stop("Cluster-robust variance estimation will not work when the data only includes a single cluster.")
  X_list <- clubSandwich:::matrix_list(X, cluster, "row")
  W_list <- clubSandwich:::weightMatrix(obj, cluster)
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, 
                 w = W_list)
  if (is.null(target)) {
    if (inverse_var) {
      Theta_list <- lapply(W_list, function(w) MASS::ginv(w))
    }
    else {
      Theta_list <- clubSandwich:::targetVariance(obj, cluster)
    }
  }
  else {
    if (!is.list(target)) {
      if (length(target) != N) {
        target <- clubSandwich:::handle_vectors(target, obj)
      }
      Theta_list <- clubSandwich:::matrix_list(target, cluster, "both")
    }
    else {
      Theta_list <- target
    }
  }

  if (type %in% c("CR2", "CR4")) {
    S <- augmented_model_matrix(obj, cluster, inverse_var, 
                                ignore_FE)
    if (is.null(S)) {
      rm(S)
      U_list <- X_list
      UW_list <- XW_list
    }
    else {
      U <- cbind(X, S)
      rm(S)
      U_list <- clubSandwich:::matrix_list(U, cluster, "row")
      UW_list <- Map(function(u, w) as.matrix(t(u) %*% 
                                                w), u = U_list, w = W_list)
    }
    UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, 
                    u = U_list)
    M_U <- matrix_power(Reduce("+", UWU_list), p = -1)
  }
  # fix this line
  # Specify the package and function name
  package_name <- "clubSandwich"
  
  # Get the function from the package
  if (type == "CR3") {
    adjustments <- do.call(pseudo.CR3, args = mget(names(formals(pseudo.CR3))))
  } else {
    adjustments <- do.call(type, args = mget(names(formals(type))))
  }

  E_list <- clubSandwich:::adjust_est_mats(type = type, est_mats = XW_list, 
                            adjustments = adjustments)
  resid <- clubSandwich:::residuals_CS(obj)
  res_list <- split(resid, cluster)
  components <- do.call(cbind, Map(function(e, r) e %*% r, 
                                   e = E_list, r = res_list))
  v_scale <- clubSandwich:::v_scale(obj)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) 
    w_scale <- 1L
  if (form == "estfun") {
    bread <- sandwich::bread(obj)
    estfun <- bread %*% components
    return(estfun * (w_scale/v_scale))
  }
  meat <- tcrossprod(components) * w_scale^2/v_scale
  if (form == "sandwich") {
    bread <- sandwich::bread(obj)
  }
  else if (form == "meat") {
    bread <- NULL
  }
  else if (is.matrix(form)) {
    bread <- form
    form <- "sandwich"
  }
  vcov <- switch(form, sandwich = bread %*% meat %*% bread/v_scale, 
                 meat = meat)
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "type") <- type
  attr(vcov, "cluster") <- cluster
  attr(vcov, "bread") <- bread
  attr(vcov, "v_scale") <- v_scale
  attr(vcov, "est_mats") <- XW_list
  attr(vcov, "adjustments") <- adjustments
  attr(vcov, "target") <- Theta_list
  attr(vcov, "inverse_var") <- inverse_var
  attr(vcov, "ignore_FE") <- ignore_FE
  class(vcov) <- c("vcovCR", "clubSandwich")
  return(vcov)
}

pseudo.CR3 <- function (X_list, XW_list) 
{
  # Multiply XW_list with X_list element-wise
  XWX_list <- Map(function(xw, x) xw %*% x, xw = XW_list, x = X_list)
  
  # Sum the matrices in XWX_list
  XWX_sum <- Reduce("+", XWX_list)
  
  # Replace chol2inv with generalized inverse using MASS::ginv
  
  M <- MASS::ginv(XWX_sum) # Use generalized inverse instead of chol2inv
  
  # Call IH_jj_list function
  IH_jj <- clubSandwich:::IH_jj_list(M, X_list, XW_list)
  
  # Apply generalized inverse instead of solve
  lapply(IH_jj, ginv)
}
