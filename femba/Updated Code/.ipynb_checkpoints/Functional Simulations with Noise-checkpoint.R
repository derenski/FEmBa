### Try increasing T
### Try increasing K 
library(stringr)
library(reshape2)
library(openxlsx)
library(mvtnorm)
library(splines2)
library(fastICA)
library(ggplot2)
library(LaplacesDemon)
library(MASS)
library(glmnet)
library(stringr)
library(dplyr)
library(gtools)
library(splines)
library(foreach)
library(doParallel)

options( warn = -1 )

literal_eval <- function(input, allowed = c("list", "c", "+", "-", "/", "*", "rep")) {
  # Create safe empty environment
  safe_env <- new.env(parent = emptyenv())
  
  # assign allowed functions
  lapply(allowed,function(f) assign(f,get(f, "package:base"),safe_env))
  
  # Evaluate input
  safe_env$expr <- parse(text = input)
  eval(substitute(expr,env = safe_env), env = safe_env)
}


desired_clusters <- 8

max_available_clusters <- detectCores()

cl <- makeCluster(min(c(max_available_clusters, desired_clusters))-3)

registerDoParallel(cl)

#library(LogConcDEAD) ## Estimation of log-concave densities with MLE approach
### Broken package


#x <- matrix(rnorm(10000),ncol=8)
#mledens   <- mlelcd(x=x)





source("Function Data Simulation Parameter Makers.R")
source("Tweedie's Formula Multivariate and Score Function.R")
source("Eigenbasis Estimator With Traditional PCA.R")

##### FOR UPLOADING
## # scp -r ./Desktop/USC/research/functional_data_tweedie/CODE/'For Uploading To CRAN'/'functions and sim code'
# derenski@msbjlr04.marshall.usc.edu:/home/derenski/femba 

## FOR DOWNLOADING
# rsync -avh derenski@msbjlr04.marshall.usc.edu:/home/derenski/femba/reports ./Desktop/USC/research/functional_data_tweedie/CODE/'For Uploading To CRAN'

nodeName <- Sys.info()["nodename"] 

sim_params <- read.xlsx("functional_simulation_parameters.xlsx", sheet="Tweedie Simulations",
                        rowNames = T)

chosen_prior <- sim_params["chosen_prior",]

distributionParameters <- literal_eval(sim_params['prior_parameters',])

full_n <- as.numeric(sim_params['N',])

p_mu <- as.numeric(sim_params['p_mu',])

p_gamma <- as.numeric(sim_params['p_g',])

iterations <- as.numeric(sim_params['iterations',])

rho_mu = as.numeric(sim_params['rho_mu',])

rho_g = as.numeric(sim_params['rho_g',])

sigma_e <- as.numeric(sim_params['sigma_e',])

SNR <- as.numeric(sim_params['snr',])

perfect_covariance <- as.logical(sim_params["perfect_covariance", ])

curveModelingMethod <- sim_params["curve modeling method", ]

Sigma_mu <- make_rho_mat(rho = rho_mu, p = p_mu)

Sigma_gamma <- make_rho_mat(rho = rho_g, p = p_gamma)

smallest_log_lambda = -6

largest_log_lambda = 6

cumu_var = .99


take_top <- as.logical(sim_params["take_top", ])

top_k_values <- ceiling(as.numeric(sim_params['top_k',])*full_n) ## Look at most extreme, where we are likely to perform best








### TO DO
### Fix the result visualization portion of the code

all_function_names <- as.vector(lsf.str())

functional_name <- all_function_names[str_detect(all_function_names, pattern="^tau_")]



set.seed(144)
###### Autocovariance structure (p^n)



### Do a Beta, or Mixture prior to get FEmBa to perform best 

######################################################## SIMULATIONS BEGIN
#########################################################




thetas <- c()

times <- seq(.01,1,.01)

mu_m <- rep(0,p_mu)

## You can add more priors if you desire, you can also specify parameters 
## in the excel file if you wish (get that working tomorrow).



es <- diag(rep(1,p_mu))

scale_mat <- diag(c(.5,7,3,1))

size_of_eigenbasis <- as.list(rep(NA, iterations))

## Define error object here

error_data <- array(NA, dim=c(6, 9, 5, iterations))

if(perfect_covariance){
  
  covariance_distance <- rep(0, length.out=iterations)
  
}else{

  covariance_distance <- seq(0,1, length.out=iterations)

}

ourPrior <- qarb(do.call(chosen_prior, c(list(n=10000), distributionParameters)))

############################################################### POINT OF NO RETURN
for (j in 1:iterations){
  
  error_array_this_iteration <- array(NA, dim=c(6, 9, 5))
  
  ## Note: you can get vectors by using mean vectors instead of scalars 
  
  simulation_parameters <- curve_generator(n_obs=full_n, Sigma_mu=Sigma_mu, Sigma_gamma=Sigma_gamma,
                                           SNR=SNR, 
                                           sigma_e=sigma_e, unmixedCoordinateDist=ourPrior,
                                           times=times)
  
  ms <- t(simulation_parameters$S_mu %*% simulation_parameters$mu)
  
  W_theta <- simulation_parameters$W_theta
  
  gs <- t(simulation_parameters$S_gamma %*% simulation_parameters$gamma)
  
  xs <- ms+gs+simulation_parameters$white_noise
  
  cov_gamma <- simulation_parameters$cov_gamma
  
  smoother <- function(D,V){
    
    smoothed <- smooth.spline(x = D, y = V)
    
    return(smoothed$y)
    
  }

  curves_for_smoothing <- as.list(data.frame(t(xs)))
  
  times_for_eigenbasis <- lapply(apply(matrix(rep(times, each=full_n), nrow=full_n), 
                                       MARGIN=1, FUN= list), FUN=unlist)
  
  smoothed_curves <- t(apply(xs, MARGIN=1, FUN=smoother, D = times))
  
  listed_smoothed_curves <- lapply(apply(smoothed_curves, MARGIN=1, FUN= list),
                                   FUN=unlist)

  
  if (curveModelingMethod == 'FPCA'){

  fpca_information <- basis_and_coefficient_estimator_pca(X = listed_smoothed_curves,
                                                  our_times = times_for_eigenbasis, 
                                                  sample_splitting = F,
                                                  cumulative_variation_described = cumu_var,
                                                  number_of_comps = 2, 
                                                  choose_comp_by_proportion = T
                                                 )

  modeling_basis <-  fpca_information$first_sample$eigenbasis_for_curves
  
  size_of_eigenbasis <- c(size_of_eigenbasis, dim(modeling_basis)[2])
  
  #### ONLY CHECK SAMPLED CURVES, CODE UNDER CONSTRUCTION!!! 
  #### (finish subsampling in "basis_estimator")
  
  correctionSample <- str_extract(fpca_information$first_sample$curve_ids_we_can_correct_with_this_split,
                                   pattern = "[0-9]+")
  
  }else{
    
    chosen_df <- dimensionalityCalculator(xs=times_for_eigenbasis,             ys=listed_smoothed_curves, pGrid=seq(5, 10, 1))   
      
      
    modeling_basis <- bSpline(times, df=chosen_df, intercept = T)
    
    correctionSample <- seq(1, full_n, 1)
    
  }
  
  correctionSample <- as.numeric(correctionSample)
  
  n_obs <- length( correctionSample)
  
  the_mus <- simulation_parameters$mu[, correctionSample]
  
  actual_curves <- ms[correctionSample, ]
  
  estimated_thetas <- thetaCalculator(yList = listed_smoothed_curves[correctionSample], 
                                basis=modeling_basis)
  
  fpca_smoothed_curves <- smoothed_curves[correctionSample, ]
  
  unsmoothed_curves <- xs[correctionSample, ]
  
  fpca_smoothed_curves <- t(modeling_basis %*% estimated_thetas)
  
  estimated_sigma_square <- mean(
    rowMeans((unsmoothed_curves-fpca_smoothed_curves)^2))
  
  newParameters <- parameterTransformerNewBasis(basis=modeling_basis, Sigma_mu = Sigma_mu, true_mu_basis = simulation_parameters$S_mu,
                                                Sigma_g=cov_gamma, true_g_basis=simulation_parameters$S_gamma,
                                                sigma_e=sigma_e)
  
  Sigma_mu_basis <- newParameters$Sigma_mu_basis

  Sigma_g_basis <- newParameters$Sigma_g_basis
   
  Sigma_e_basis <- newParameters$Sigma_e_basis

   estimated_conditional_sigma <- Sigma_g_basis+Sigma_e_basis
    
  perturbed_covariance <- covariance_sampler(estimated_conditional_sigma, 
                                               desired_distance=0)

  final_ests_trans <- tweedieCorrectionNonJointICA(X=estimated_thetas, 
                                     Sigma_gamma=round(perturbed_covariance, 3),  
                                     functionalBasis=modeling_basis,  
                                     transformation="decorrelate",
                                     gridStep = .001, 
                                     lambdaGrid=10^seq(-6, 6, 1),
                                     numberOfFolds = 5, 
                                     numberOfKnotsScoreFunction = 8)
  
  print('First FEmBA')
  
  # truehist(the_mus-final_ests_trans[[1]])

  thetas_tweedie_trans <- final_ests_trans$tweedieEstimates
  
  #   conditional_sigma = Sigma_g
  #    x = estimated_thetas
  #      transformation = T
  #    lambda = lambd
  #    step = .001
  #      data_basis = bases_from_data  
  
  final_ests_notrans <- tweedieCorrectionNonJointICA(X=estimated_thetas, 
                                     Sigma_gamma=round(perturbed_covariance, 3),  
                                     functionalBasis=modeling_basis,  
                                     transformation="none",
                                     gridStep = .001, 
                                     lambdaGrid=10^seq(-6, 6, 1),
                                     numberOfFolds = 5, 
                                     numberOfKnotsScoreFunction = 8)
    
    
    
  
  print('Second FEmBA')
  
  thetas_tweedie_ICA <- tweedieCorrectionWithICA(X=estimated_thetas, varianceCovariance=round(perturbed_covariance, 3), 
                                                 functionalBasis=modeling_basis,
                                                 gridStep = .001, 
                                                 lambdaGrid=10^seq(-6,6,1),
                                                 numberOfFolds = 5, 
                                                 numberOfKnotsScoreFunction = 8,
                                                 maxIterations=100)
  
  ## plot((thetas_tweedie_ICA-estimated_thetas)[2,] ~ estimated_thetas[2,])
  

  thetas_tweedie_notrans <- final_ests_notrans$tweedieEstimates
  
  #plot((thetas_tweedie_notrans-estimated_thetas)[2,] ~ estimated_thetas[2,])
  
  thetas_js_trans <- james_stein(X=estimated_thetas, 
          relevant_Sigma_g=round(perturbed_covariance, 3), 
          transform = T)
  
  
  thetas_js_notrans <- james_stein(X=estimated_thetas, 
                                    relevant_Sigma_g=round(perturbed_covariance, 3), 
                                    transform = F)
  ### Used when first choosing prior to get marginal. Now we choose marginal to get prior
  thetas_oracle <- oracle_estimator(the_thetas=estimated_thetas, 
                                      Sigma_mu=Sigma_mu_basis, Sigma_g=Sigma_g_basis,
                                      the_mus=least_squares(modeling_basis, y=simulation_parameters$S_mu %*% the_mus))
  
 # Sigma_gamma_forOracle <- covariance_basis_transformer(current_covariance=cov_gamma,
  #                                                      old_basis=simulation_parameters$S_gamma,
  #                                                      new_basis= modeling_basis)
  
  
  #thetas_oracle <- tweediesFormulaOracle(X=estimated_thetas, W=W_theta, U =simulation_parameters$U,
  #                                       Sigma_gamma=Sigma_gamma_forOracle, 
  #                                       functionalBasis=modeling_basis)
  
  ## plot((thetas_oracle-simulation_parameters$theta)[1,] ~ simulation_parameters$theta[1,])

  tweedie_curves_trans <- t(modeling_basis %*% thetas_tweedie_trans)
    
  tweedie_curves_notrans <- t(modeling_basis %*% thetas_tweedie_notrans)
  
  tweedie_curves_ICA <- t(modeling_basis %*% thetas_tweedie_ICA)
    
  js_curves_trans <- t(modeling_basis %*% thetas_js_trans)
    
  js_curves_notrans <- t(modeling_basis %*% thetas_js_notrans)

  oracle_curves_notrans <- t(modeling_basis %*% thetas_oracle)

  all_curves <- list(actual_curves, tweedie_curves_trans, tweedie_curves_notrans, tweedie_curves_ICA,
                     js_curves_trans, js_curves_notrans,
                     oracle_curves_notrans,
                     fpca_smoothed_curves, xs)
  
  names(all_curves) <- c('actual', "tweedie_trans", "tweedie_notrans", "tweedie_ICA","js_trans",
                         "js_notrans", "oracle", "smoothed", "unsmoothed")
  
  for (func_slice in 1:length(functional_name)){
  
    ordering_functional <- get(functional_name[func_slice])
    
    all_vals <- lapply(all_curves, FUN=function(x) apply(x, MARGIN=1, 
                                                         FUN=ordering_functional))
    
    smoothed_vals <- all_vals[["smoothed"]]
    
    vals_of_interest <- lapply(all_vals, FUN=function(x) x[order(smoothed_vals, 
      decreasing = take_top)[1:top_k_values]])
    
    error_data_val <- sapply(vals_of_interest, FUN=bias_variance_mse, 
                             y=vals_of_interest[["actual"]])
    
    error_array_this_iteration[, , func_slice] <- error_data_val
  
  }
  
  ##### L2 difference, make sure this works
  
  curve_l2_norms <- lapply(all_curves, FUN = function(x) apply(x, MARGIN=1, 
                        FUN=function(x) mean(x^2)))
  
  l2_smoothed <- curve_l2_norms[["smoothed"]]
  
  curves_of_interest <- lapply(all_curves, FUN=function(x) x[order(l2_smoothed, 
        decreasing = T)[1:top_k_values], ]) ## Look for largest initial norm
  
  l2_dist_data <- sapply(curves_of_interest, FUN=bias_variance_mse, 
                         y=curves_of_interest[["actual"]])
  
  error_array_this_iteration[, , dim(error_array_this_iteration)[3]] <- l2_dist_data
  
  dimnames(error_array_this_iteration) <- list(rownames(l2_dist_data),
                                               colnames(l2_dist_data),
                    c(functional_name, "l2 distance"))
  
  if(j==1){
    
    dimnames(error_data)[1:3] <- dimnames(error_array_this_iteration)
    
    dimnames(error_data)[[4]] <- paste("iteration_", seq(1, iterations, 1),
                                     sep="")
    
  }
  
  error_data[, , , j] <- error_array_this_iteration
  
  # dim1: error summary; dim2: correction method; dim3: functional; dim4: iteration
  
  print(j)
  
}


#apply(error_data, MARGIN=c(1,2,3), FUN=mean)
























report_dir <- paste("../reports",nodeName, sep="/")


if (!exists(report_dir)){
  
  dir.create(report_dir, recursive = T)
  
}

date_dir <- paste(report_dir, Sys.Date(), sep="/")

if (!dir.exists(date_dir)){
  
  dir.create(date_dir)
  
}

output_dir <- paste(date_dir, paste("sim",length(list.files(date_dir))+1, chosen_prior, "prior", "SNR", SNR, 
                                    c("perturbation", "performance")[perfect_covariance+1], "analysis", 
                                    c("bottom", "top")[take_top+1], 100*as.numeric(sim_params['top_k',]), 
                                    "percent" , sep="_"), sep="/")

plot_dir <- paste(output_dir, "plots", sep='/')

table_dir <- paste(output_dir, "tables", sep='/')

if (!dir.exists(output_dir)){
  
  dir.create(output_dir)
  
  dir.create(plot_dir)
  
}






## Get score function example
score_f_js <- thetas_js_trans-estimated_thetas

score_f_femba <- thetas_tweedie_trans-estimated_thetas


(ggplot(NULL, aes(x=estimated_thetas[1,], y=score_f_js[1,])) + geom_point() +
    geom_point(aes(x=estimated_thetas[1,], y=score_f_femba[1,]), color='red'))














### Check method performance on the extremes as well
apply(error_data, MARGIN=c(1, 2, 3), FUN=mean, na.rm=T)



# dim1: error summary; dim2: correction method; dim3: functional; dim4: iteration

# Capitalizing function
capitalizer <- function(word){
  
  parsed_word <- unlist(strsplit(word, split=""))
  
  parsed_word[1] <- toupper(parsed_word[1])
  
  return(paste(parsed_word, sep="", collapse=""))
  
}

metric_names <- dimnames(error_data)[[1]]

performance_metric_of_interest <- metric_names[3]

if (perfect_covariance){
  
  dir.create(table_dir)
  
  performance_plot_dir <- paste(plot_dir, "performance_plots", sep="/")
  
  dir.create(performance_plot_dir)

  for (f in c(functional_name, "l2 distance")){
    
    functional_plot_dir <- paste(performance_plot_dir, f, sep="/")
    
    dir.create(functional_plot_dir)
    
    functional_for_plotting <- f
    
    metric_plots <- as.list(rep(NA, dim(error_data)[1]))
    
    names(metric_plots) <- dimnames(error_data)[[1]]
    
    for (metric in dimnames(error_data)[[1]]){
    
      performance_metric_of_interest <- metric
      
      this_error_data_plotting <- t(error_data[performance_metric_of_interest,-1, 
                                               functional_for_plotting, ])
      
      df <- melt(this_error_data_plotting ,  id.vars = rownames, variable.name = 'Method')
      
      names(df) <- c("iteration", "Method", performance_metric_of_interest)
      
      methods <- c(expression("Unsmooth"),expression("Smooth"), expression(JS[NT]), 
                   expression(JS[T]), expression(FEmBa[NT]), expression(FEmBa[T]),expression(FEmBa["ICA"]),
                   expression("Oracle"))
      
      names(methods) <- c("unsmoothed", "smoothed", "js_notrans", "js_trans", 
                          "tweedie_notrans", "tweedie_trans", "tweedie_ICA", "oracle")
      
      current_plot <- (ggplot(df, aes_string(x="Method", y=performance_metric_of_interest, 
                             fill="Method")) + geom_boxplot()+theme_bw(base_size = 20) +
          scale_x_discrete(labels=methods[df$Method]) 
        + scale_fill_manual(
          labels=methods[df$Method],
          values=gg_color_hue(length(methods)))
        +ylab(capitalizer(performance_metric_of_interest)) +
          ggtitle(paste("Distrbution of", capitalizer(performance_metric_of_interest),
                        "by Method")))
      
      ggsave(filename=paste(functional_plot_dir, metric, sep="/"), plot=current_plot,
             width=10, height=10, units='in', dpi=300, device="pdf")
      
    
    }
    
  }
    
  
  
  pairs <- list(c('bias', "variance"), c('mse', "var_mse"), c("mape", "var_mape"))
  
  for (p in pairs){
  
    method_ordering <- c("oracle", "tweedie_trans", "tweedie_notrans", "js_trans",
                         "js_notrans", "smoothed", "unsmoothed")
    
    functional_ordering <- c("l2 distance", "tau_max", "tau_var", "tau_beta", "tau_mu")
    
    functional_performance_tables <- apply(error_data, MARGIN=c(1, 2, 3), FUN=mean, na.rm=T)
    
    this_table_data <- functional_performance_tables[p[1], -1,]
    
    se_this_table_data <- functional_performance_tables[p[2], -1,]
    
    final_mse_data <- round(this_table_data, 4)[method_ordering, functional_ordering]
    
    if (!(p[1]=='bias')){
    
    final_se_data <- round(sqrt(se_this_table_data), 4)[method_ordering, functional_ordering]
    
    }else{
      
      final_se_data <- round(se_this_table_data, 4)[method_ordering, functional_ordering]
      
    }
    
    combined_product <- array(NA, dim=dim(final_se_data))
    
    dimnames(combined_product) <- list(rownames(final_se_data), colnames(final_se_data))
    
    combined_product[] <- paste(final_mse_data, ' (', final_se_data, ')', sep="")
    
    rownames(combined_product) <- c("Oracle", "FEmBa with back-transformation",
                                    "FEmBa with no back-transformation", 
                                    "James-Stein with back-transformation",
                                    "James-Stein with no back-transformation",
                                    "Smoothed Data",
                                    "Unsmoothed Data")
    
    combined_product <- cbind(rep("", dim(combined_product)[1]), rownames(combined_product), 
                              combined_product)
    
    rownames(combined_product) <- NULL
    
    colnames(combined_product)[1:2] <- c("Prior", "Method")
    
    colnames(combined_product)[3:length(colnames(combined_product))] <- c("$L^{2}$ Norm", 
          "\\\\LARGE\\{$\\\\taumax$\\}", "\\\\LARGE\\{$\\\\tauvar$\\}",
          "\\\\LARGE\\{$\\\\taubeta$\\}", "\\\\LARGE\\{$\\\\taumu$\\}")
    
    combined_product[3:4, 1] <- c(capitalizer(chosen_prior), paste("SNR=", round(SNR, 1)
                                                                   , sep=""))
    
    latexed_names <- paste(colnames(combined_product), collapse = " & ")
    
    table_pure <- knitr::kable(combined_product, format = 'latex')
    
    table_character <- as.character(table_pure)
    
    
    correct_table_style <- "\\{ll|ccccc\\}"
    
    fixed_table <- str_replace(table_character, pattern="(?<=\\\\begin\\{tabular\\})\\{.*\\}", 
                replacement =correct_table_style)
    
    fixed_table <- str_replace(fixed_table, pattern = "Prior.*",
                               replace=latexed_names)
    
    fixed_table <- str_replace(fixed_table, pattern="\\\\hline(?=\n & Oracle)", replace="\\\\hline\n\\\\hline")
    
    fileConn<-file(paste(table_dir, paste(p[1], " and ", p[2],".txt", sep=""), sep="/"))
    writeLines(fixed_table, fileConn)
    close(fileConn)
  
  }

}else{
  
  perturbation_plot_dir <- paste(plot_dir, "perturbation_plots", sep="/")
  
  dir.create(perturbation_plot_dir)
  
  options(scipen = 999)
  
  iterations_and_distances <- data.frame(cbind(dimnames(error_data)[[4]], covariance_distance),
                                         stringsAsFactors = F)
  
  iterations_and_distances$covariance_distance <- as.numeric(
    iterations_and_distances$covariance_distance)
  
  names(iterations_and_distances)[1] <- "Iteration"
  
  
  for (f in c(functional_name, "l2 distance")){
    
    functional_plot_dir <- paste(perturbation_plot_dir, f, sep="/")
    
    dir.create(functional_plot_dir)
    
    for (metric_name in dimnames(error_data)[[1]]){
      
      perturbation_data_for_plot <- t(error_data[metric_name, -1, f, ])
      
      df <- reshape2::melt(perturbation_data_for_plot ,  
                           id.vars = rownames, variable.name = 'Method')
      
      df[,1] <- as.character(df[,1])
      df[,2] <- as.character(df[,2])
      
      names(df) <- c("Iteration", "Method", "metric_name")
      
      df_merged <- df %>% inner_join(iterations_and_distances)
      
      final_df <- df_merged %>% filter(Method %in% c("tweedie_trans", "js_trans", "smoothed"))
      
      plot_perturb <- (ggplot(final_df, 
                              aes_string(x = "covariance_distance", y = "metric_name", 
                                         col="Method")) +
                        geom_smooth(method='loess', se=F) +theme_bw(base_size = 12)+
                       xlab('Distance') + ylab(capitalizer(metric_name)) 
                       + scale_color_grey(name = "Method",
                                          labels = c(expression(FEmBa[T]),
                                                     expression(JS[T]),
                                                     "Smoothed"))
                       + theme(legend.position='none'))
      
      ggsave(filename=paste(functional_plot_dir, metric_name, sep="/"), plot=plot_perturb,
             width=10, height=10, units='in', dpi=300, device="pdf")
      
    }
    
  }
  
}

write.csv(sim_params, file=paste(output_dir, "simulation_parameters.csv", sep="/"))

#########################################################################################


  