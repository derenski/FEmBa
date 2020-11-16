library(splines2)
library(tidyr)
library(caret)

### Coefficient vector estimator 






coeff_vec <- function(y, observed_basis){ ## Estimator for basis coefficients using ols
  
  final_coeffs <- (solve((t(observed_basis) %*% observed_basis)) %*% 
                     t(observed_basis)) %*% y ##### Estimate of coefficient vector
  
  return(final_coeffs)
}

thetaCalculator <- function(yList, basis){
  
  thetas <- array(  unlist(lapply(yList, coeff_vec, observed_basis=basis)), dim=c(dim(basis)[2], length(yList)))
  
  return(thetas)
  
}

parameterTransformerNewBasis <- function(basis, Sigma_mu, true_mu_basis, Sigma_g, true_g_basis, sigma_e){
  
  Sigma_mu_basis <- covariance_basis_transformer(current_covariance=Sigma_mu, 
                                                old_basis=true_mu_basis,
                                                new_basis=basis)
  
  Sigma_mu_basis <- .5*(Sigma_mu_basis+t(Sigma_mu_basis))
  
  Sigma_g_basis <- covariance_basis_transformer(current_covariance=Sigma_g, 
                                                old_basis=true_g_basis,
                                                new_basis=basis)
  
  Sigma_g_basis <- .5*(Sigma_g_basis+t(Sigma_g_basis))
  
  Sigma_e_basis <- covariance_basis_transformer(
    current_covariance=sigma_e*diag(rep(1, dim(basis)[1])), 
    old_basis=diag(rep(1, dim(basis)[1])),
    new_basis=basis)
  
  Sigma_e_basis <- .5*(Sigma_e_basis+t(Sigma_e_basis))
  
  Sigma_theta_basis <- Sigma_g_basis + Sigma_e_basis
  
  outputList <- list(Sigma_mu_basis, Sigma_g_basis, Sigma_theta_basis, Sigma_e_basis)
  
  names(outputList) <- c("Sigma_mu_basis", "Sigma_g_basis", "Sigma_theta_basis", "Sigma_e_basis")
  
  return(outputList)
  
}









#### Function used within the 'basis_estimator_pca' function
#### Responsible for the actual PCA and loadings estimation 

fpca_estimator_trad <- function(curves_for_estimating_eigenfunctions,
                                curves_to_be_modeled, times_for_curves_to_be_modeled, 
                                variation_described = .95, n_comps = 6,
                                choose_comps_with_proportion = T, 
                                number_of_iterations = 1000,
                                we_want_a_constraint = F,
                                transformation_we_are_using = NULL){

  fpca_results <- princomp(x=curves_for_estimating_eigenfunctions)
  
  cum_var <- cumsum(fpca_results$sdev^2)/sum(fpca_results$sdev^2)
  
  if (choose_comps_with_proportion == T){
  
  chosen_size = min(which(cum_var > variation_described))
  
  } else if (choose_comps_with_proportion == F){
    
    chosen_size = n_comps
    
  }
  
  basis <- fpca_results$loadings[,1:chosen_size]

  Si_s <- list()
  
  ### It's already done the scaling, you don't need the Ds
  ### Try Just (S^T)X
  
  mats_to_save <- list()
  
  unique_ids_for_curves_we_can_correct <- row.names(curves_to_be_modeled)
  
  theta_iks <- matrix(NA, nrow = chosen_size, 
                      ncol = length(unique_ids_for_curves_we_can_correct))
  
  min_times <- as.numeric(unlist(lapply(times_for_curves_to_be_modeled, FUN = min, na.rm = T)))
  
  max_times <- as.numeric(unlist(lapply(times_for_curves_to_be_modeled, FUN = max, na.rm = T)))
  
  begin_time <- min(min_times)
  
  end_time <- max(max_times)
  
  full_grid <- times_for_curves_to_be_modeled[[1]]
  
  domain <- times_for_curves_to_be_modeled[[1]]
  
  basis_matrix <- matrix(NA, nrow = length(domain), 
                           ncol = chosen_size)
    
  for (i in 1:chosen_size){
      
      mod1 <- splinefun(x = full_grid, y = basis[,i],
                        method = 'natural')
      
      basis_matrix[,i] <- mod1(domain)
      
  }
  
  S_theta <- basis_matrix
  
  final_stuff <- list(S_theta, 
                      row.names((curves_to_be_modeled)))
  
  names(final_stuff) <- c( "eigenbasis_for_curves",
                           "curve_ids_we_can_correct_with_this_split")
  
  return(final_stuff)
  
}


#########################################################
################ Below is the function used to perform PCA
########### On the observed functional data.
# Only use this code for simulations, or if your curves
### are densely observed, with no gaps. 


basis_and_coefficient_estimator_pca <- function(X, sample_splitting = F,
            our_times = list(), cumulative_variation_described = .9,
            number_of_comps = 6, choose_comp_by_proportion = T){
  
  n_obs <- length(X)
  
  ids <- paste0("curve_",1:n_obs)
  
  curve_data_final <- matrix(unlist(X), nrow=n_obs, byrow=T)
  
  row.names(curve_data_final) <- ids
  
  
  if (sample_splitting){
    
    curve_samp <- sample(1:n_obs, size = n_obs, replace = F)
    
    samp_1 <- curve_samp[1:(.5*n_obs)]
    
    samp_1 <- samp_1[order(samp_1, decreasing = F)]
    
    samp_2 <- curve_samp[(.5*n_obs+1):n_obs]
    
    samp_2  <- samp_2 [order(samp_2, decreasing = F)]
    
    split_1 <- curve_data_final[samp_1,]
    
    split_2 <- curve_data_final[samp_2 ,]
    
    times_for_split_1 <- our_times[samp_1]
    
    times_for_split_2 <- our_times[samp_2]
    
    
  }else{
    
    ids_for_estimating_eigenfunctions <- ids
    
    split_1 <- curve_data_final
    
    split_2 <- curve_data_final
    
    times_for_split_1 <- our_times
    
    times_for_split_2 <- our_times
    
  } 

  if (sample_splitting){
    
    results_for_split_1 <- fpca_estimator_trad(curves_for_estimating_eigenfunctions = split_2,
                                               curves_to_be_modeled = split_1,
                                               times_for_curves_to_be_modeled = times_for_split_1,
                                               variation_described = cumulative_variation_described,
                                               n_comps = number_of_comps,
                                               choose_comps_with_proportion = choose_comp_by_proportion)
    
    results_for_split_2 <- fpca_estimator_trad(curves_for_estimating_eigenfunctions = split_1,
                                               curves_to_be_modeled = split_2,
                                               times_for_curves_to_be_modeled = times_for_split_2,
                                               variation_described = cumulative_variation_described,
                                               n_comps = number_of_comps,
                                               transformation_we_are_using = transformation_for_constraints)
    
    all_needed_info <- list(results_for_split_1, results_for_split_2)
    
    names(all_needed_info) <- c('first_sample','second_sample')
    
    return(all_needed_info)
    
  }else{
    
    results_for_split_1 <- fpca_estimator_trad(curves_for_estimating_eigenfunctions = split_2,
                                               curves_to_be_modeled = split_1,
                                               times_for_curves_to_be_modeled = times_for_split_1,
                                               variation_described = cumulative_variation_described,
                                               n_comps = number_of_comps,
                                               choose_comps_with_proportion = choose_comp_by_proportion)
    
    all_needed_info <- list(results_for_split_1)
    
    names(all_needed_info) <- c('first_sample')
    
    return(all_needed_info)
    
  }
  
}

######################## Dimensionality estimation with pre-chosen basis

### xs: list of domain values, ys: list of range values
dimensionalityCalculator <- function(xs, ys, pGrid=seq(5, 10, 1)){
  
  bicMatrix <- array(NA, dim=c(length(xs), length(pGrid)))
  
  dimnames(bicMatrix) <- c(list(seq(1,length(xs), 1)), list(pGrid))
  
  for (curve in 1:length(xs)){
  
    this_x <- xs[[curve]]
    
    this_y <- ys[[curve]]
      
    train_x <- sample(this_x, size=floor(.7*length(this_x)), replace=F)
      
    train_points <- (this_x %in% train_x)
      
    test_points <- !train_points
  
    special_bases <- lapply(pGrid, 
                              FUN = function(i) bSpline(x = this_x,
                                                        degree = 3, intercept = T, 
                                                        df = i))

    model_fitter <- function(design_matrix, y, train_points){

        model_vector <- coeff_vec(y[train_points], observed_basis=design_matrix[train_points,])
        
        predictions <- design_matrix[test_points,] %*% model_vector
        
        the_mse <- mean((predictions-y[test_points])^2)

        return(the_mse)

      }

    complexity_values <- unlist(lapply(special_bases, FUN = model_fitter, y = this_y, train_points=train_points))

    names(complexity_values) <- pGrid
                            
    bicMatrix[curve,] <- complexity_values
  
  }
  
  best_values <- as.numeric(dimnames(bicMatrix)[[2]][apply(bicMatrix, MARGIN=1, FUN=which.min)])
                            
  return(floor(mean(best_values)))
   
}
  
  


