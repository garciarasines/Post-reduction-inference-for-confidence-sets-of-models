rm(list = ls())
pkgs = c("glmnet", "mvtnorm", "future.apply", "pracma")
installed = pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed])
invisible(lapply(pkgs, library, character.only = TRUE))

source("~/.../Gaussian_functions.R")

two_stage = function(n, p, theta, ks, sigma, rho, max_size_model = 5, gamma = 0.6, alpha = 0.05, B = 500){
  
  # Correct model 
  E_null = which(theta != 0)        
  
  # 1 if encompassing set includes true model
  sure_screening_COX = numeric(B)             
  sure_screening_split_COX = numeric(B)
  sure_screening_LASSO = numeric(B)             
  sure_screening_split_LASSO = numeric(B)
  
  # 1 if true model belongs to the confidence set
  cov_COX_Q1 = numeric(B)                     
  cov_COX_Q2 = numeric(B)
  cov_COX_R = numeric(B)  
  cov_COX_LRT = numeric(B)    
  cov_COX_LRTsplit = numeric(B)    
  cov_LASSO_Q1 = numeric(B)
  cov_LASSO_Q2 = numeric(B)
  cov_LASSO_R = numeric(B)
  cov_LASSO_LRT = numeric(B)
  cov_LASSO_LRTsplit = numeric(B)
  
  # size of confidence set
  size_COX_Q1 = numeric(B)                    
  size_COX_Q2 = numeric(B)
  size_COX_R = numeric(B)
  size_COX_LRT = numeric(B)    
  size_COX_LRTsplit = numeric(B)
  size_LASSO_Q1 = numeric(B)                    
  size_LASSO_Q2 = numeric(B)
  size_LASSO_R = numeric(B)
  size_LASSO_LRT = numeric(B)    
  size_LASSO_LRTsplit = numeric(B)
  
  pb = txtProgressBar(min = 0, max = 100, style = 3, width = 50, char = "=") 
  
  # Simulation
  set.seed(12345)
  for(b in 1:B){
    
    # Data
    X = design(n, p, rho)
    y = y_sample(X, theta, sigma)
    sigma_hat = sigma_estimator(y, X, gamma)
    
    # Encompassing models
    E_hat_COX = encompassing("Cox", y, X)
    sure_screening_COX[b] = all(E_null %in% E_hat_COX)
    E_hat_split_COX = encompassing("Cox", y[1:(0.6*n)], X[1:(0.6*n), ])
    sure_screening_split_COX[b] = all(E_null %in% E_hat_split_COX)
    
    E_hat_LASSO = encompassing("lasso", y, X)
    sure_screening_LASSO[b] = all(E_null %in% E_hat_LASSO)
    E_hat_split_LASSO = encompassing("lasso", y[1:(0.6*n)], X[1:(0.6*n), ])
    sure_screening_split_LASSO[b] = all(E_null %in% E_hat_split_LASSO)
    
    # 1st confidence set using Q 
    cs_COX_Q1 = confset("Q", y, X, sigma_hat, E_hat_COX, max_size_model, alpha, ks[1], gamma)
    cov_COX_Q1[b] = any(sapply(cs_COX_Q1, function(x) identical(x, E_null)))
    size_COX_Q1[b] = length(cs_COX_Q1)  
    
    cs_LASSO_Q1 = confset("Q", y, X, sigma_hat, E_hat_LASSO, max_size_model, alpha, ks[1], gamma)
    cov_LASSO_Q1[b] = any(sapply(cs_LASSO_Q1, function(x) identical(x, E_null)))
    size_LASSO_Q1[b] = length(cs_LASSO_Q1) 
    
    # 2nd confidence set using Q 
    cs_COX_Q2 = confset("Q", y, X, sigma_hat, E_hat_COX, max_size_model, alpha, ks[2], gamma)
    cov_COX_Q2[b] = any(sapply(cs_COX_Q2, function(x) identical(x, E_null)))
    size_COX_Q2[b] = length(cs_COX_Q2)  
    
    cs_LASSO_Q2 = confset("Q", y, X, sigma_hat, E_hat_LASSO, max_size_model, alpha, ks[2], gamma)
    cov_LASSO_Q2[b] = any(sapply(cs_LASSO_Q2, function(x) identical(x, E_null)))
    size_LASSO_Q2[b] = length(cs_LASSO_Q2)  
    
    # Confidence set using R
    cs_COX_R = confset("R", y, X, sigma_hat, E_hat_COX, max_size_model, alpha, k = NULL, gamma)
    cov_COX_R[b] = any(sapply(cs_COX_R, function(x) identical(x, E_null)))
    size_COX_R[b] = length(cs_COX_R)  
    
    cs_LASSO_R = confset("R", y, X, sigma_hat, E_hat_LASSO, max_size_model, alpha, k = NULL, gamma)
    cov_LASSO_R[b] = any(sapply(cs_LASSO_R, function(x) identical(x, E_null)))
    size_LASSO_R[b] = length(cs_LASSO_R)  
    
    # Confidence set using LRT
    cs_COX_LRT = confset("LRT", y, X, sigma_hat, E_hat_COX, max_size_model, alpha, k = NULL, gamma)
    cov_COX_LRT[b] = any(sapply(cs_COX_LRT, function(x) identical(x, E_null)))
    size_COX_LRT[b] = length(cs_COX_LRT)  
    
    cs_LASSO_LRT = confset("LRT", y, X, sigma_hat, E_hat_LASSO, max_size_model, alpha, k = NULL, gamma)
    cov_LASSO_LRT[b] = any(sapply(cs_LASSO_LRT, function(x) identical(x, E_null)))
    size_LASSO_LRT[b] = length(cs_LASSO_LRT) 
    
    # Confidence set using LRT with sample splitting
    cs_COX_LRTsplit = confset("LRT", y[-(1:(0.6*n))], X[-(1:(0.6*n)), ], sigma_hat, E_hat_split_COX, max_size_model, alpha, k = NULL, gamma)
    cov_COX_LRTsplit[b] = any(sapply(cs_COX_LRTsplit, function(x) identical(x, E_null)))
    size_COX_LRTsplit[b] = length(cs_COX_LRTsplit) 
    
    cs_LASSO_LRTsplit = confset("LRT", y[-(1:(0.6*n))], X[-(1:(0.6*n)), ], sigma_hat, E_hat_split_LASSO, max_size_model, alpha, k = NULL, gamma)
    cov_LASSO_LRTsplit[b] = any(sapply(cs_LASSO_LRTsplit, function(x) identical(x, E_null)))
    size_LASSO_LRTsplit[b] = length(cs_LASSO_LRTsplit)  
    
    setTxtProgressBar(pb, 100*b/B)
    
  }
  
  results = list()
  
  # Cox results
  results[[1]] = rbind(c(mean(sure_screening_COX), mean(sure_screening_split_COX)), c(sd(sure_screening_COX)/sqrt(B), sd(sure_screening_split_COX)/sqrt(B)))
  results[[2]] = rbind(c(mean(cov_COX_Q1), mean(cov_COX_Q2), mean(cov_COX_R), mean(cov_COX_LRT), mean(cov_COX_LRTsplit)),
                       c(sd(cov_COX_Q1)/sqrt(B), sd(cov_COX_Q2)/sqrt(B), sd(cov_COX_R)/sqrt(B), sd(cov_COX_LRT)/sqrt(B), sd(cov_COX_LRTsplit)/sqrt(B)))
  results[[3]] = rbind(c(mean(size_COX_Q1), mean(size_COX_Q2), mean(size_COX_R), mean(size_COX_LRT), mean(size_COX_LRTsplit)),
                       c(sd(size_COX_Q1)/sqrt(B), sd(size_COX_Q2)/sqrt(B), sd(size_COX_R)/sqrt(B), sd(size_COX_LRT)/sqrt(B), sd(size_COX_LRTsplit)/sqrt(B)))
  
  # Lasso results
  results[[4]] = rbind(c(mean(sure_screening_LASSO), mean(sure_screening_split_LASSO)), c(sd(sure_screening_LASSO)/sqrt(B), sd(sure_screening_split_LASSO)/sqrt(B)))
  results[[5]] = rbind(c(mean(cov_LASSO_Q1), mean(cov_LASSO_Q2), mean(cov_LASSO_R), mean(cov_LASSO_LRT), mean(cov_LASSO_LRTsplit)),
                       c(sd(cov_LASSO_Q1)/sqrt(B), sd(cov_LASSO_Q2)/sqrt(B), sd(cov_LASSO_R)/sqrt(B), sd(cov_LASSO_LRT)/sqrt(B), sd(cov_LASSO_LRTsplit)/sqrt(B)))
  results[[6]] = rbind(c(mean(size_LASSO_Q1), mean(size_LASSO_Q2), mean(size_LASSO_R), mean(size_LASSO_LRT), mean(size_LASSO_LRTsplit)),
                       c(sd(size_LASSO_Q1)/sqrt(B), sd(size_LASSO_Q2)/sqrt(B), sd(size_LASSO_R)/sqrt(B), sd(size_LASSO_LRT)/sqrt(B), sd(size_LASSO_LRTsplit)/sqrt(B)))
  
  return(results)
  
}

# Simulation parameters
p = 400
ks = c(2, 8)
sigma = 1
ns = c(100, 150)
ts = c(0.5, 1)
rhos = c(0.1, 0.5)

# Example: table 1
n = ns[2]
t = ts[2]
rho = rhos[2]
theta = c(rep(t, 3), rep(0, p - 3))

two_stage(n, p, theta, ks, sigma, rho) 
