gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

make_rho_mat <- function(rho,p){
  
  the_vec <- matrix(NA, nrow = p, ncol =  p)
  
  for (i in 1:(p)){
    
    for (j in 1:p){
      
      the_vec[i,j] <- rho^abs(i-j)
      
    }
  }
  
  return(.5*(the_vec+t(the_vec)))
  
}


dens_calculator <- function(X,MU, precision_matrix){
  
  # return(dmvnorm(x = X, mean = MU, sigma = condit_cov, log = F))
  
  X_mat_form <- matrix(X, ncol =  1)
  
  quad_form_thing <- -.5*(t(X_mat_form-MU) %*% precision_matrix %*% (X_mat_form-MU))
  
  dens_val <- (det((1/(2*pi)) * precision_matrix)^(1/2))*exp(quad_form_thing)
  
  return(dens_val)
  
  
}

dens_calculator_vectorized <- function(X,MU, precision_matrix){
  
  # return(dmvnorm(x = X, mean = MU, sigma = condit_cov, log = F))
  
  X_mat_form <- matrix(rep(X, dim(MU)[2]), nrow = dim(MU)[1])
  
  quad_form_thing <- -.5*(t(X_mat_form-MU) %*% precision_matrix %*% (X_mat_form-MU))
  
  dens_mat <- (det((1/(2*pi)) * precision_matrix)^(1/2))*exp(quad_form_thing)
  
  return(diag(dens_mat))
  
  
}


dens_calculator_vectorized_over_X <- function(X,MU, condit_cov){
  
  dens_vals <- dmvn(X, mu = MU, Sigma = condit_cov)
  
  return(dens_vals)
  
}







gradient_calculator <- function(x_given_mu, mu, 
                                inverted_conditional_sigma, dens_val){
  
  dens_coeff <- dens_val
  
  return(dens_coeff*inverted_conditional_sigma %*% (mu-x_given_mu) )
  
}

gradient_calculator_vectorized <- function(x_given_mu, mu, 
                                           inverted_conditional_sigma, dens_coeffs){
  
  dens_coeff <- matrix(dens_coeffs, ncol =  1)
  
  x_for_mult <- matrix(rep(x_given_mu, dim(mu)[2]), nrow = dim(mu)[1])
  
  first_part <- inverted_conditional_sigma %*% (mu-rep(x_given_mu) )
  
  gradient_estimates <- sapply(1:length(dens_coeff), FUN = function(i) dens_coeff[i] * first_part[,i] )
  
  return(gradient_estimates)
  
}

point_mass_score_estimator <- function(ps,point_masses, precision_matrix,x){
  
  density_values <- dens_calculator_vectorized(X = x, MU = point_masses, 
                                               precision_matrix = precision_matrix)
  
  density_values <- matrix(density_values, ncol =  1)
  
  gradients <- gradient_calculator_vectorized(x_given_mu = x, mu = point_masses,
                                              inverted_conditional_sigma = precision_matrix,
                                              dens_coeffs = density_values)
  
  numerator <- gradients %*% ps
  
  denominator <- as.numeric((ps) %*% density_values)
  
  return(numerator/denominator)
  
}





mixture_score_estimator <- function(ps,point_masses, covariance_matrices,x){
  
  density_values <- as.numeric(lapply(1:dim(point_masses)[2], FUN = function(i) dmvn(x, 
                    mu = point_masses[,i], Sigma = covariance_matrices[[i]]) ))
  
  density_values <- as.matrix(density_values, ncol =  1)
  
  precision_mats <- lapply(covariance_matrices, FUN = solve)
  
  gradients_list_form <-lapply(1:dim(point_masses)[2], 
                               FUN = function(i) gradient_calculator(x_given_mu = x, mu = point_masses[,i],
                                                                     inverted_conditional_sigma =precision_mats[[i]],
                                                                     dens_val = as.numeric(density_values[i,])) )
  
  gradients <- matrix(NA,nrow = lengths(gradients_list_form)[1], 
                      ncol =  length(gradients_list_form))
  
  for (grad in 1:length(gradients_list_form)){
    
    gradients[,grad] <- gradients_list_form[[grad]]
    
  }
  
  
  numerator <- gradients %*% ps
  
  denominator <- as.numeric(t(ps) %*% density_values)
  
  return(numerator/denominator)
  
}

# Curve generator
qarb <- function(emp_dist){
  
  quantile_func <- function(p){
  
    initial_dist <- quantile(emp_dist, p, na.rm=TRUE)
    
    initial_dist[is.infinite(initial_dist)] <- 0
    
    return(initial_dist)
  
  }
  
  return(quantile_func)
  
}


norta <- function(n_obs, corr_mat, distribution=qnorm){
  
  Zs <- rmvn(n=n_obs, mu=rep(0, dim(corr_mat)[1]), Sigma=corr_mat)
  
  phi_Zs <- pnorm(Zs)
  
  phi_Zs[phi_Zs %in% c(0,1)] <- .5
  
  desired_variables <- apply(phi_Zs, MARGIN=1, FUN=distribution)
   
  return(desired_variables)
  
}




### New data generator, infer prior from marginal
curve_generator <- function(n_obs, misspessMat, Sigma_mu, Sigma_gamma, SNR, 
                            sigma_e, unmixedCoordinateDist,
                            times, xSatisfyICA=FALSE){
  
  ICAmodelDataMaker <- function(n, p, quantile_func, misspessMat, covMat){
    
    ### Generate unmixed data with NORTA
    XUnmixed <- norta(n, corr_mat = misspessMat, distribution = quantile_func)
      
    XUnmixed <- matrixToPower(cov(t(XUnmixed)), -.5) %*% XUnmixed
    
    ### True U in so(p)
      ### Try generatingU from different distributions rexp(dim(XUnmixed )[1]^2, rate=10)
    UTrue <- gramSchmidt(array(runif(dim(XUnmixed )[1]^2, min=-1, max=1), 
                               dim=rep(dim(XUnmixed )[1], 2)))$Q 
    
    ### Covariance matrix to the -1/2 power
    covMatMinusHalf <- matrixToPower(covMat, -.5)
    
    ### True unmixing matrix
    trueW <- UTrue %*% covMatMinusHalf
    
    ### The observed, mixed data
    X <- solve(trueW) %*% XUnmixed
    
    finalOutput <- list(X=X, W=trueW)
    
    return(finalOutput)
    
  }
  
  p_mu <- dim(Sigma_mu)[1]
  
  p_gamma <- dim(Sigma_gamma)[1]
  
  S_mu <- bSpline(times, degree = min(p_mu,3), df = p_mu, intercept = T)
  
  S_gamma <- bSpline(times, degree = min(p_gamma,3), df = p_gamma, intercept = T)
    
  if (xSatisfyICA){
      
      
      if (dim(Sigma_mu)[1] >= dim(Sigma_gamma)[1]){
      
          Sigma_gamma <- (solve(t(S_gamma) %*% S_gamma) %*% t(S_gamma) %*% S_mu %*% 
                          Sigma_mu %*% t(S_mu) %*%  S_gamma %*% solve(t(S_gamma) %*% S_gamma) )


          Sigma_gamma <- .5*(Sigma_gamma+t(Sigma_gamma))
          
          }else{
          
          Sigma_mu <- (solve(t(S_mu) %*% S_mu) %*% t(S_mu) %*% S_gamma %*% 
                          Sigma_gamma %*% t(S_gamma) %*%  S_mu %*% solve(t(S_mu) %*% S_mu) )


          Sigma_mu <- .5*(Sigma_mu+t(Sigma_mu))
          
          
      }
      
  }
  
  raw_snr <- (norm(Sigma_mu,type = 'F')
              /norm(Sigma_gamma, type = 'F'))
  
  muModelData <- ICAmodelDataMaker(n=n_obs, p=p_mu, quantile_func=unmixedCoordinateDist,                                            misspessMat=misspessMat, covMat=Sigma_mu)
  
  mus <- muModelData$X
  
  Wtheta <- muModelData$W
  
  U <- muModelData$W %*% matrixToPower(Sigma_mu, .5)
  
  gammas <- (raw_snr/SNR)*norta(n_obs=n_obs, corr_mat=Sigma_gamma, distribution=qnorm)
  
  Xs <- t(S_mu %*% mus + S_gamma %*% gammas)

  e_noise <- rmvnorm(n_obs, mean = rep(0,length(times)),
                     sigma  = sigma_e*diag(length(times)))  
  
  final_output <- list(mu=mus, U=U, S_mu=S_mu, Xs=Xs, gamma=gammas,  S_gamma=S_gamma, white_noise=e_noise)
  
  return(final_output)
  
}

## Covariace basis transformer

covariance_basis_transformer <- function(current_covariance, old_basis,
                                         new_basis){
  
  transform_matrix <- solve(t(new_basis) %*% new_basis) %*% t(new_basis)
 
  estimated_conditional_sigma_inner_term <- (old_basis %*% current_covariance %*% 
                                               t(old_basis))
  
  return(transform_matrix %*% 
      estimated_conditional_sigma_inner_term %*% t(transform_matrix))
   
}







##### Function for perturbing covariance matrix

covariance_sampler <- function(estimated_conditional_sigma, desired_distance=0){

  p <- dim(estimated_conditional_sigma)[1]
  
  for_perturbing_matrix <- rWishart(1, df= .5*p*(p+1), 
                                    Sigma=estimated_conditional_sigma) 
  
  noise_for_covariance <- ((apply(for_perturbing_matrix, MARGIN=1, FUN=rowMeans)/(.5*p*(p+1)))
                           -estimated_conditional_sigma)
  
  ## Makes noise exactly symmetric
  
  noise_for_covariance <- .5*(noise_for_covariance+t(noise_for_covariance))
  
  initial_distance <- norm(solve(estimated_conditional_sigma) %*% 
                             noise_for_covariance, '2')
  
  noise_for_covariance <- desired_distance*noise_for_covariance/initial_distance
  
  return(estimated_conditional_sigma+noise_for_covariance)

}


### Various functionals of interest 

tau_mu <- function(x){
  
  return(mean(x, na.rm=T))
  
  
}

tau_max <- function(x){
  
  return(max(x, na.rm=T))
  
  
}

tau_var <- function(x){
  
  return(max(x, na.rm=T))
  
  
}

### Various competitor methods

james_stein <- function(X, relevant_Sigma_g, transform=F){
  
  sigma <- cov(t(X))
  
  p <- dim(sigma)[1]
  
  svdecomp <- svd(sigma)
  
  sigma_minus_half <- (svdecomp$u) %*% diag(svdecomp$d^-.5) %*% t(svdecomp$v)
  
  sigma_one_half <- (svdecomp$v) %*% diag(svdecomp$d^.5) %*% t(svdecomp$u)
  
  trans_X <- sigma_minus_half %*% X
    
  if (transform){
      
    js_ests <- ((diag(rep(1, p)) - sigma_minus_half %*% 
                   relevant_Sigma_g %*% sigma_minus_half) %*% trans_X
                + sigma_minus_half %*% relevant_Sigma_g %*% sigma_minus_half %*% 
                  t(sapply(rowMeans(trans_X), FUN=rep, dim(X)[2])))
    
    js_ests <- sigma_one_half %*% js_ests
    
  }else{
    
    the_scale_est <- solve(sigma)
    
    js_ests <- ((diag(rep(1, p)) - relevant_Sigma_g %*% the_scale_est) %*% X
                    + relevant_Sigma_g %*% the_scale_est %*% 
                      t(sapply(rowMeans(X), FUN=rep, dim(X)[2])))
    
    }
  
  return(js_ests)
  
}


least_squares <- function(X, y){
  
  return(solve(t(X) %*% X) %*% t(X) %*% y)
  
}

oracle_estimator <- function(the_thetas, 
        Sigma_mu, Sigma_g, the_mus){
    
    conditional_density_estimates <- t(apply(the_mus, MARGIN = 2, 
                                             FUN = function(x) dens_calculator_vectorized_over_X(
                                               X= t(the_thetas), condit_cov = Sigma_g, MU = x)))
    
    marginal_density_estimates <- apply(conditional_density_estimates, MARGIN =2,
                                        FUN = mean)
    
    conditional_precision <- solve(Sigma_g)
    
    conditional_gradients <- lapply(1:dim(the_thetas)[2], FUN = function(i) 
      gradient_calculator_vectorized(x_given_mu = the_thetas[,i], mu = the_mus,
                                     inverted_conditional_sigma = conditional_precision,
                                     dens_coeffs= conditional_density_estimates[,i]))
    
    gradient_estimates <- lapply(conditional_gradients, FUN = rowMeans)
    
    gradient_estimates_final <- matrix(unlist(gradient_estimates),
                                       nrow = dim(Sigma_g)[1])
    
    score_function_estimates <- sapply(1:dim(gradient_estimates_final)[2],
                                       FUN = function(i) gradient_estimates_final[,i]/marginal_density_estimates[i]  )
    
    oracle_thetas <- the_thetas + (Sigma_g) %*% score_function_estimates
    
    return(oracle_thetas)
  
}





##### Functionals for curves

tau_mu <- function(curve){
  
  return(mean(curve, na.rm=T))

  }

tau_max <- function(curve){
  
  return(max(curve, na.rm=T))

  }

tau_var <- function(curve){

  return(var(curve, na.rm=T))

  }

tau_beta <- function(curve, beta_t = function(t) t){
  
  times <- seq(.01, 1, .01)
  
  curve_matrix_form <- matrix(curve, nrow=1, ncol=length(curve))
  
  beta_matrix_form <- as.matrix(beta_t(times), nrow=length(times), ncol=1)
  
  return((curve_matrix_form %*% beta_matrix_form)/length(times))
  
}

bias_variance_mse <- function(y_hat, y){
  
  y_hat <- as.numeric(y_hat)
    
  y <- as.numeric(y)  
    
  bias <- mean(y_hat-y)
  
  variance <- sd(y_hat-y)^2
  
  mse <- mean((y_hat-y)^2)
  
  var_mse <- (sd((y_hat-y)^2)^2)/length(y)
  
  mape <- mean(abs((y_hat-y)/y), na.rm=T)
  
  var_mape <- (sd(abs((y_hat-y)/y), na.rm=T)^2)/length(y)
  
  final_output <- c(bias, variance, mse, var_mse, mape, var_mape)
  
  names(final_output) <- c('bias', "variance", "mse", "var_mse", "mape",
                           "var_mape")
  
  return(final_output)
  
}


