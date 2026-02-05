# Data generation
design = function(n, p, rho, Sigma_str = "toeplitz"){
  if (Sigma_str == "toeplitz"){
    Sigma = outer(1:p, 1:p, function(i, j) rho^abs(i-j))     
  } else if (Sigma_str == "block"){
    Sigma = diag(p)
    Sigma[1:10, 1:10] = (1 - rho)*diag(10) + rho*matrix(1, 10, 10)
  } 
  return(rmvnorm(n, mean = rep(0, p), sigma = Sigma))
}
y_sample = function(X, theta, sigma){
  n = nrow(X)
  return(X%*%theta + sigma*rnorm(n))
}

# Randomization coefficients
coefs = function(k, sigma_hat){
  as = numeric(k-1)
  bs = numeric(k-1)
  as[1] = sigma_hat*sqrt(k-1)
  bs[1] = sigma_hat^2/as[1]
  if (k > 2){
    for (i in 2:(k-1)){
      as[i] = sqrt((k-1)*sigma_hat^2-sum(bs[1:(i-1)]^2))
      bs[i] = (sigma_hat^2+sum(bs[1:(i-1)]^2))/as[i]
    }
  }
  mat = matrix(0, nrow = k, ncol = k)
  mat[1, ] = rep(1, k)
  for(i in 2:k){
    mat[i, i-1] = as[i-1]
    mat[i, i:k] = -bs[i-1]
  }
  return(mat)
}

# Construction of encompassing model
Cox_square = function(y, X, alpha_cox){
  
  M = 1:ncol(X)
  
  # Square size 
  ss = ceiling(sqrt(p))
  total_cells = ss^2
  
  # Randomly arrange indices into the cube; missing cells filled with NAs
  idx = sample(M)
  idx = c(idx, rep(NA, total_cells - p))  
  square = matrix(idx, nrow = ss, ncol = ss)
  
  selected_indices = c()
  
  # Function to run regression and return variables significant at level alpha
  significant_vars = function(y, X_sub, vars, alpha) {
    df = data.frame(y, X_sub)
    fit = lm(y ~ . - 1, data = df)
    pvals = summary(fit)$coefficients[ , 4]  
    sig_idx = which(pvals < alpha)
    if (length(sig_idx) == 0){
      return(integer(0))
    } else{
    return(vars[sig_idx])
    }
  }
  
  # Iterate over rows and columns
  for (i in 1:ss) {
    # Row
    vars_row = na.omit(square[i, ])
    if (length(vars_row) >= 1) {
      X_sub = X[ , vars_row, drop = FALSE]
      sig = significant_vars(y, X_sub, vars_row, alpha_cox)
      selected_indices = c(selected_indices, sig)
    }
    # Column
    vars_col = na.omit(square[, i])
    if (length(vars_col) >= 1) {
      X_sub = X[, vars_col, drop = FALSE]
      sig = significant_vars(y, X_sub, vars_col, alpha_cox)
      selected_indices = c(selected_indices, sig)
    }
  }
  
  S_hat = sort(unique(selected_indices))
  return(S_hat)
}
Cox = function(y, X, max_size = 15){

    alpha_cox = 0.051
    E_hat = 1:(max_size + 1)
    
    while(length(E_hat) > max_size){
      alpha_cox = alpha_cox - 0.001
      E_hat = Cox_square(y, X, alpha_cox)
    }
  
  return(E_hat)
}
lasso = function(y, X, max_size = 15){
  fit = glmnet(X, y, intercept = FALSE)
  df = fit$df
  lambdas = fit$lambda
  lambda_selected = min(lambdas[df <= max_size])
  coef_vec = as.vector(coef(fit, s = lambda_selected))
  E_hat = which(coef_vec[-1] != 0)
  return(E_hat)
}
encompassing = function(type, y, X){
  
  if (type == "Cox"){
    E_hat = Cox(y, X)
  } else if (type == "lasso"){
    E_hat = lasso(y, X)
  }
  return(E_hat)
}

# Sd estimation
sigma_estimator = function(y, X, gamma, type = "mrcv"){
  
  n = length(y)
  n2 = floor(gamma*n)
  if (n2 %% 2 != 0) n2 = n2 + 1
  
  y_new = y[1:n2]
  X_new = X[1:n2, ]
  
  y1 = y_new[1:(n2/2)]
  y2 = y_new[-(1:(n2/2))]
  X1 = X_new[1:(n2/2), ]
  X2 = X_new[-(1:(n2/2)), ]
  
  # First estimator
  S = encompassing("lasso", y2, X2)
  df1 = n2/2 - length(S)
  X1_S = X1[ , S]
  M1 = diag(1, n2/2) - X1_S%*%solve(t(X1_S)%*%X1_S)%*%t(X1_S)
  sigma2_hat_1 = sum((M1%*%y1)^2)/(n2/2-length(S))
  
  # Second estimator
  S = encompassing("lasso", y1, X1)
  df2 = n2/2 - length(S)
  X2_S = X2[ , S]
  M2 = diag(1, n2/2) - X2_S%*%solve(t(X2_S)%*%X2_S)%*%t(X2_S)
  sigma2_hat_2 = sum((M2%*%y2)^2)/(n2/2 - length(S))
  
  # Average of estimators
  if (type == "rcv"){
    sigma_hat = sqrt(0.5*(sigma2_hat_1 + sigma2_hat_2))
  } else if (type == "mrcv"){
    sigma_hat = sqrt((df1*sigma2_hat_1 + df2*sigma2_hat_2)/(df1 + df2))
  }
  
  return(sigma_hat)
  
}

# Submodel test
pval = function(type, y, y_splits, X, sigma_hat, E_hat, E_null, k, gamma){
  
  n = nrow(X)
  d = length(E_null)
  
  X_null = X[ , E_null]
  M = diag(1, n) - X_null%*%solve(t(X_null)%*%X_null)%*%t(X_null)
  
  if (type == "Q"){
    
    U = eigen(M)$vectors[ , 1:(n-d)]
    
    # Matrix of "q" copies, stored in columns
    qs_unnorm = t(U)%*%y_splits
    qs = qs_unnorm%*%diag(1/sqrt(diag(t(qs_unnorm)%*%qs_unnorm)))

    inner_prods = t(qs)%*%qs
    p_val = pnorm((sqrt(2*(n-d))/k)*sum(inner_prods[upper.tri(inner_prods)]), lower.tail = FALSE)
    
  } else if (type == "R"){
    
    R2 = sum((M%*%y)^2) 
    
    p_val = pchisq(R2/sigma_hat^2, df = n - d, lower.tail = FALSE)
    
  } else if (type == "LRT"){
    
    fit_enc = lm(y ~ X[ , E_hat] - 1)
    fit_null = lm(y ~ X[ , E_null] - 1)
    lrt = anova(fit_enc, fit_null)
    p_val = lrt$`Pr(>F)`[2]    
    
  } 
  
  return(p_val)
  
}

# Computation of confidence set
confset = function(type, y, X, sigma_hat, E_hat, max_size_model, alpha, k, gamma){
  
  # Split of y for Q-test: same split is used for all submodels' tests
  y_splits = NULL
  if (type ==  "Q"){
    y_splits = cbind(y, matrix(rnorm(n*(k-1)), nrow = n, ncol = k-1))%*%coefs(k, sigma_hat)
  }
  
  # Test every submodel of size =< max_size
  confset = list()
  for (size in 1:max_size_model){
    subsets = combn(E_hat, size, simplify = FALSE)
    for (E_null in subsets) {
      p_val = pval(type, y, y_splits, X, sigma_hat, E_hat, E_null, k, gamma)
      if (p_val > alpha) confset[[length(confset) + 1]] = E_null
    }
  }
  
  return(confset)

}








