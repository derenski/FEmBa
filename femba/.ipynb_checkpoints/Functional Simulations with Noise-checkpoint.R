### Try increasing T
### Try increasing K 
library(stringr)
library(reshape2)
library(openxlsx)

### TO DO
### Fix the result visualization portion of the code

date_dir <- paste(".", Sys.Date(), sep="/")

if (!dir.exists(date_dir)){
  
  dir.create(date_dir)
  
}

output_dir <- paste(date_dir, paste("sim_",length(list.files(date_dir))+1, sep=""),
                    sep="/")

plot_dir <- paste(output_dir, "plots", sep='/')

table_dir <- paste(output_dir, "tables", sep='/')

if (!dir.exists(output_dir)){
  
  dir.create(output_dir)
  
  dir.create(plot_dir)
  
  dir.create(table_dir)
  
  
}

source("Function Data Simulation Parameter Makers.R")
source("Tweedie's Formula Multivariate and Score Function.R")
source("Eigenbasis Estimator With Traditional PCA.R")

all_function_names <- as.vector(lsf.str())

functional_name <- all_function_names[str_detect(all_function_names, pattern="^tau_")]



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
set.seed(144)
###### Autocovariance structure (p^n)



### Do a Beta, or Mixture prior to get FEmBa to perform best 

######################################################## SIMULATIONS BEGIN
#########################################################


sim_params <- read.xlsx("functional_simulation_parameters.xlsx", sheet=1,
                        rowNames = T)

chosen_prior <- sim_params["chosen_prior",]

full_n <- as.numeric(sim_params['N',])

p_mu <- as.numeric(sim_params['p_mu',])

p_gamma <- as.numeric(sim_params['p_g',])

iterations <- as.numeric(sim_params['iterations',])

rho_mu = as.numeric(sim_params['rho_mu',])

rho_g = as.numeric(sim_params['rho_g',])

sigma_e <- as.numeric(sim_params['sigma_e',])

SNR <- as.numeric(sim_params['snr',])

Sigma_mu <- make_rho_mat(rho = rho_mu, p = p_mu)

Sigma_g <- make_rho_mat(rho = rho_g, p = p_gamma)
  
raw_snr <- norm(Sigma_mu,type = 'F')/norm(Sigma_g, type = 'F')

Sigma_g <- (raw_snr/SNR)*Sigma_g

smallest_log_lambda = -6

largest_log_lambda = 6

cumu_var = .99

top_k_values <- ceiling(as.numeric(sim_params['top_k',])*full_n) ## Look at most extreme, where we are likely to perform best

thetas <- c()

times <- seq(.01,1,.01)

mu_m <- rep(0,p_mu)

## You can add more priors if you desire, you can also specify parameters 
## in the excel file if you wish (get that working tomorrow).


if (chosen_prior == "mixture"){

   our_points <- cbind(sample(c(1,2), size = p_mu, replace = T), 
                 sample(c(-1,-2), size = p_mu, replace = T))

   probs <- rep(1/dim(our_points)[2], dim(our_points)[2])
   
   sgn <- c(-1, 1)[rbern(1000, .5)+1]
   
   our_dist <- rnorm(10000, mean=sgn*2) 
   
   
}else if (chosen_prior == "gaussian"){
  
  our_dist <- rnorm(10000)
  
}else if (chosen_prior=="gamma"){
  
  our_dist <- rgamma(10000, shape=1, scale=5)
  
}else{
  
  stop('Specify a valid distribution')
  
}

qarb <- function(p, emp_dist){

  quantile(emp_dist, p)
  
}

es <- diag(rep(1,p_mu))

scale_mat <- diag(c(.5,7,3,1))

lambdas <- rep(NA, iterations)

cv_errs_trans <- as.list(rep(NA, iterations))

cv_errs_notrans <- as.list(rep(NA, iterations))

size_of_eigenbasis <- as.list(rep(NA, iterations))

## Define error object here

error_data <- array(NA, dim=c(6, 8, 5, iterations))

covariance_distances <- seq(0,1, length.out=iterations)

############################################################### POINT OF NO RETURN
for (j in 1:iterations){
  
  error_array_this_iteration <- array(NA, dim=c(6, 8, 5))
  
  n <- full_n
  
  ## Note: you can get vectors by using mean vectors instead of scalars 
  
  simulation_parameters <- curve_generator(n_obs=1000, Sigma_mu=Sigma_mu, Sigma_g=Sigma_g, 
                                           sigma_e=sigma_e, prior=qarb, emp_dist=our_dist,
                                           times=seq(.01, 1, .01))
  
  ms <- t(simulation_parameters$m_basis %*% simulation_parameters$mu)
  
  gs <- t(simulation_parameters$g_basis %*% simulation_parameters$gamma)
  
  xs <- ms+gs+simulation_parameters$white_noise
  
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


  fpca_information <- basis_and_coefficient_estimator_pca(X = listed_smoothed_curves,
                                                  our_times = times_for_eigenbasis, 
                                                  sample_splitting = F,
                                                  cumulative_variation_described = cumu_var,
                                                  number_of_comps = 4, 
                                                  choose_comp_by_proportion = T
                                                 )
  
  all_the_bases <- fpca_information$first_sample$eigenbases_for_curves
  
  all_the_transforms <- fpca_information$first_sample$transforms_for_gamma

  fpca_basis_for_integral <-  fpca_information$first_sample$basis_for_integral_calc
  
  size_of_eigenbasis <- c(size_of_eigenbasis, dim(fpca_basis_for_integral)[2])

  fpca_basis <-  fpca_information$first_sample$eigenbases_for_curves[[1]]
  
  # fpca_basis <- s_basis
  
  #### ONLY CHECK SAMPLED CURVES, CODE UNDER CONSTRUCTION!!! 
  #### (finish subsampling in "basis_estimator")
  
  curve_ids_we_can_correct_with_this_split <- str_extract(fpca_information$first_sample$curve_ids_we_can_correct_with_this_split,
                                   pattern = "[0-9]+")
  
  curve_ids_we_can_correct_with_this_split <- as.numeric(curve_ids_we_can_correct_with_this_split )
  
  n_obs <- length( curve_ids_we_can_correct_with_this_split)
  
  the_mus <- simulation_parameters$mu[, curve_ids_we_can_correct_with_this_split]
  
  actual_curves <- ms[curve_ids_we_can_correct_with_this_split, ]
  
  the_thetas <- fpca_information$first_sample$all_thetas[,
              curve_ids_we_can_correct_with_this_split ]
  
  fpca_smoothed_curves <- smoothed_curves[curve_ids_we_can_correct_with_this_split, ]
  
  unsmoothed_curves <- xs[curve_ids_we_can_correct_with_this_split, ]

  
  fpca_smoothed_curves <- t(fpca_information$first_sample$eigenbases_for_curves[[1]] %*%
                             the_thetas)
  
  estimated_sigma_square <- mean(
    rowMeans((unsmoothed_curves-fpca_smoothed_curves)^2))
  
   Sigma_mu_fpca <- covariance_basis_transformer(current_covariance=Sigma_mu, 
                      old_basis=simulation_parameters$m_basis,
                      new_basis=fpca_basis)
   
   
   Sigma_g_fpca <- covariance_basis_transformer(current_covariance=Sigma_g, 
                                                 old_basis=simulation_parameters$g_basis,
                                                 new_basis=fpca_basis)
   
   Sigma_e_fpca <- covariance_basis_transformer(
     current_covariance=sigma_e*diag(rep(1, length(times))), 
    old_basis=diag(rep(1, length(times))),
      new_basis=fpca_basis)
   
   estimated_conditional_sigma <- Sigma_g_fpca+Sigma_e_fpca
  
  # covariance_distances[j]  
    
  perturbed_covariance <- covariance_sampler(estimated_conditional_sigma, 
                                               desired_distance=0)

  final_ests_trans <- tweedie_correction(X = the_thetas, 
                      variance_covariance = round(perturbed_covariance, 3),
                      transformation = T, step = .001, ICA = F,
                      functional_basis = fpca_basis_for_integral,
                      lambda_grid=10^seq(-6, 6, 1),
                      number_of_folds = 5, number_of_knots_score_function = 8)
  
  print('First FEmBA')
  
  # truehist(the_mus-final_ests_trans_1[[1]])
  
  
  cv_errs_trans[[j]] <- final_ests_trans[[3]]
  
  lambdas[j] <- final_ests_trans[[2]]
  
  thetas_tweedie_trans <- final_ests_trans[[1]]
  
  #   conditional_sigma = Sigma_g
  #    x = the_thetas
  #      transformation = T
  #    lambda = lambd
  #    step = .001
  #      data_basis = bases_from_data  
  
  final_ests_notrans <- tweedie_correction(X = the_thetas, 
                        variance_covariance = round(perturbed_covariance, 3),
                        transformation = F, step = .001, ICA = F, 
                        functional_basis = fpca_basis_for_integral,
                        lambda_grid=10^seq(-6, 6, 1),
                        number_of_folds = 5, number_of_knots_score_function = 8)
  
  print('Second FEmBA')
  
  cv_errs_notrans[[j]] <- final_ests_notrans[[3]]
  
  thetas_tweedie_notrans <- final_ests_notrans[[1]]
  
  thetas_js_trans <- james_stein(X=the_thetas, 
          relevant_Sigma_g=round(perturbed_covariance, 3), 
          transform = T)
  
  
  thetas_js_notrans <- james_stein(X=the_thetas, 
                                    relevant_Sigma_g=round(perturbed_covariance, 3), 
                                    transform = F)
  
  thetas_oracle <- oracle_estimator(the_thetas=the_thetas, 
                                      Sigma_mu=Sigma_mu_fpca, Sigma_g=Sigma_g_fpca,
                                      the_mus=least_squares(fpca_basis, y=simulation_parameters$m_basis %*% the_mus))

  tweedie_curves_trans <- t(fpca_basis %*% thetas_tweedie_trans)
    
  tweedie_curves_notrans <- t(fpca_basis %*% thetas_tweedie_notrans)
    
  js_curves_trans <- t(fpca_basis %*% thetas_js_trans)
    
  js_curves_notrans <- t(fpca_basis %*% thetas_js_notrans)

  oracle_curves_notrans <- t(fpca_basis %*% thetas_oracle)

  all_curves <- list(actual_curves, tweedie_curves_trans, tweedie_curves_notrans,
                     js_curves_trans, js_curves_notrans,
                     oracle_curves_notrans,
                     fpca_smoothed_curves, xs)
  
  names(all_curves) <- c('actual', "tweedie_trans", "tweedie_notrans", "js_trans",
                         "js_notrans", "oracle", "smoothed", "unsmoothed")
  
  for (func_slice in 1:length(functional_name)){
  
    ordering_functional <- get(functional_name[func_slice])
    
    all_vals <- lapply(all_curves, FUN=function(x) apply(x, MARGIN=1, 
                                                         FUN=ordering_functional))
    
    smoothed_vals <- all_vals[["smoothed"]]
    
    vals_of_interest <- lapply(all_vals, FUN=function(x) x[order(smoothed_vals, 
      decreasing = T)[1:top_k_values]])
    
    error_data_val <- sapply(vals_of_interest, FUN=bias_variance_mse, 
                             y=vals_of_interest[["actual"]])
    
    error_array_this_iteration[, , func_slice] <- error_data_val
  
  }
  
  ##### L2 difference
  
  l2_dist_data <- sapply(all_curves, FUN=bias_variance_mse, 
                         y=all_curves[["actual"]])
  
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


## Get score function example
score_f <- thetas_tweedie_trans-the_thetas


plot(score_f[1,] ~ the_thetas[1,])

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

for (f in c(functional_name, "l2 distance")){
  
  functional_plot_dir <- paste(plot_dir, f, sep="/")
  
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
                 expression(JS[T]), expression(FEmBa[NT]), expression(FEmBa[T]),
                 expression("Oracle"))
    
    names(methods) <- c("unsmoothed", "smoothed", "js_notrans", "js_trans", 
                        "tweedie_notrans", "tweedie_trans", "oracle")
    
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


write.csv(sim_params, file=paste(output_dir, "simulation_parameters.csv", sep="/"))



if (FALSE){


##### You can fix this tomorrow. 


average_risk <- Reduce('+', cv_errs_trans)/length(cv_errs_trans)

average_risk <- data.frame(average_risk, stringsAsFactors = F)

names(average_risk) <- c('Lambda','Average.Risk')

average_risk$is_min <- c('Not Min','Min')[as.numeric(as.factor(average_risk$Average.Risk == 
                                                       
                                                                 min(average_risk$Average.Risk)))]

lambdas_trans <- c()

for (h in 1:length(cv_errs_trans)){
  
  lambdas_trans[h] <- cv_errs_trans[[h]][which.min(cv_errs_trans[[h]][,2]),1]
  }

trans_tab <- table(lambdas_trans)

mode_trans <- as.numeric(names(trans_tab)[which.max(trans_tab)])

### + scale_y_log10()\n(Log Scale)

(ggplot(average_risk, aes(x = log(Lambda)/log(10), y = Average.Risk+1)) + geom_line()

  +geom_point() +scale_x_continuous(breaks = seq(-12,6,1)) +theme_classic(base_size = 12)
  
  +xlab('Base 10 log of Lambda') + ylab('Average Cross Validated Risk')
  
  +ggtitle('Average 10-fold Cross Validated Risk, against Lambda
  
           for Transformed FEmBa')
  
  +scale_color_manual(name = 'Minimum', labels = c('The Minimum','Not the Minimum'),
                      
                      values = c('Red','Blue')) + geom_vline(xintercept = log(mode_trans)/log(10),
                                                             
                                                             col =  'red')
  )

average_risk <- Reduce('+', cv_errs_notrans)/length(cv_errs_notrans)

average_risk <- data.frame(average_risk, stringsAsFactors = F)

names(average_risk) <- c('Lambda','Average.Risk')

average_risk$is_min <- c('Not Min','Min')[as.numeric(as.factor(average_risk$Average.Risk == 
                                                                 
                                                               min(average_risk$Average.Risk)))]










lambdas_notrans <- c()

for (h in 1:length(cv_errs_trans)){
  
  lambdas_notrans[h] <- cv_errs_notrans[[h]][which.min(cv_errs_notrans[[h]][,2]),1]
}

notrans_tab <- table(lambdas_notrans)

mode_notrans <- as.numeric(names(notrans_tab)[which.max(notrans_tab)])

### + scale_y_log10() \n(Log Scale)

(ggplot(average_risk, aes(x = log(Lambda)/log(10), y = Average.Risk)) + geom_line()
  
  +geom_point() +scale_x_continuous(breaks = seq(-12,6,1)) +theme_classic(base_size = 12)
  
  +xlab('Base 10 log of Lambda') + ylab('Average Cross Validated Risk')
  
  +ggtitle('Average 10-fold Cross Validated Risk, against Lambda
           
           for Untransformed FEmBa')
  
  +scale_color_manual(name = 'Minimum', labels = c('The Minimum','Not the Minimum'),
                      
                      values = c('Red','Blue'))+ geom_vline(xintercept = log(mode_notrans)/log(10),
                                                            
                                                            col =  'red'))






#########################################################################################

options(scipen = 999)

if (metric == 'l2_norm'){

  
  
plot_perturb <- (ggplot(performance_data_for_plot, 
    aes(x = covariance_distance, y = log(l2_distance), col=method)) 
  + geom_smooth(method='loess', se=F) +theme_classic(base_size = 12)
  +xlab('Distance') + ylab('Log of average' ~ L^{2} ~ 'norm') 
  + scale_color_grey(name = "Method",
                       labels = c("Oracle", expression(FEmBa[T]),
                                  expression(FEmBa[NT]), expression(JS[T]),
                                  expression(JS[NT]), "Smooth"))
  + theme(legend.position='none'))

#  +scale_color_manual(labels = c("T999", "T888"))

}else if (metric == 'max'){

plot_perturb <- (ggplot(performance_data_for_plot, aes(x = covariance_distance, y = log(l2_distance), col =  method)) 
                 + geom_smooth(method='loess', se=F) +theme_classic(base_size = 12)
                 +xlab('Distance') + ylab('Log MSE') 
                 + ggtitle('Log MSE for Max, Across Methods')
                 + scale_color_manual(name = "Method",
                   labels = c("Oracle", expression(FEmBa[T]),
                                                 expression(FEmBa[NT]), expression(JS[T]),
                                                 expression(JS[NT]), "Smooth"),
                                      values = gg_color_hue(6))
                 + theme(legend.position='none')
)

}else if (metric == 'variance'){

plot_perturb <- (ggplot(performance_data_for_plot, aes(x = covariance_distance, y = log(l2_distance), col =  method)) 
                 + geom_smooth(method='loess', se=F) +theme_classic(base_size = 12)
                 +xlab('Distance') + ylab('Log MSE') 
                 + ggtitle('Log MSE for Variance, Across Methods')
                 + scale_color_manual(name = "Method",
                   labels = c("Oracle", expression(FEmBa[T]),
                                                 expression(FEmBa[NT]), expression(JS[T]),
                                                 expression(JS[NT]), "Smooth"),
                                      values = gg_color_hue(6))
                 + theme(legend.position='none')
)

}else if (metric == 'regression'){

plot_perturb <- (ggplot(performance_data_for_plot, aes(x = covariance_distance, y = log(l2_distance), col =  method)) 
                 + geom_smooth(method='loess', se=F) +theme_classic(base_size = 12)
                 +xlab('Distance') + ylab('Log MSE') 
                 + ggtitle('Log MSE for Regression, Across Methods')
                 + scale_color_manual(name = "Method",
                   labels = c("Oracle", expression(FEmBa[T]),
                                                 expression(FEmBa[NT]), expression(JS[T]),
                                                 expression(JS[NT]), "Smooth"),
                                      values = gg_color_hue(6))
                 + theme(legend.position='none')
)

}else if (metric == 'average'){

plot_perturb <- (ggplot(performance_data_for_plot, aes(x = covariance_distance, y = log(l2_distance), col =  method)) 
                 + geom_smooth(method='loess', se=F) +theme_classic(base_size = 12)
                 +xlab('Distance') + ylab('Log MSE') 
                 + ggtitle('Log MSE for Average Value, Across Methods')
                 + scale_color_manual(name = "Method",
                   labels = c("Oracle", expression(FEmBa[T]),
                                                 expression(FEmBa[NT]), expression(JS[T]),
                                                 expression(JS[NT]), "Smooth"),
                                      values = gg_color_hue(6))
                 + theme(legend.position='none')
)

}

plot_name <- paste0("perturbation_plot_", chosen_prior , "_prior_", metric)

ggsave(paste0("./simulation_plots/", plot_name, ".pdf"), plot=plot_perturb,
       dpi=400, width=5, height=3, units='in')





  }
  