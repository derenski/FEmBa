library(mixtools)
library(stringr) ### Loading Packages
library(knitr)
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
library(plyr)
library(dplyr)
library(gtools)
library(splines)
library(foreach)
library(doParallel)
library(clue)
library(abind)
library(pracma)
library(data.table)
library(sn)
library(xtable)


options(dplyr.summarise.inform = FALSE)

max_available_clusters <- detectCores()-1 ### Creating a cluster for parallel programming
  
desired_clusters <- 4
  
cl <- makeCluster(min(c(max_available_clusters, desired_clusters)))

registerDoParallel(cl)

tableFixer <- function(badLatex){ ### Fixes tables that have a bad latex representation
    
    fixed_table <- as.character(badLatex)
    
    fixed_table <- str_replace_all(badLatex, pattern = "textbackslash\\{\\}",
                               replace='')
      
      
    fixed_table <- str_replace_all(fixed_table, pattern = '\\\\(?=(\\{|\\}))', replace='')  
    
    return(fixed_table)
}




sciFormatter <- function(number, digits){ ### Formats numbers in scientific notation, when appropriate

        
    neededNumber <- formatC(number, format = "e", digits = digits)
        
 #   if (str_detect(neededNumber, pattern="e(\\+|-)00$")){
            
  #      neededNumber <- str_replace( neededNumber, pattern="e(\\+|-)00$", replace="")  
            
            
  #  }
        
    return(neededNumber)
        
}

averagedVectorl2Norm <- function(X,Y){ ### Calculates L2-norms between N by T matrices
                                       ### Used to calculate L2-Norms of Curves. 

    return(rowMeans((X-Y)^2))
    
}

distWithPermutAndSgnChange <- function(X, Y){
  
  indexes <- seq(1, dim(X)[1], 1)
  
  allDists <- array(NA, dim=c(length(indexes), length(indexes), 2))
  
  ### Form tensor of all possible distances between the rows of Y and the rows of X (row-permuted and signed)
  ### Dim 1: row of Y; Dim 2: row of X; Dim3: Sign of the row of X (1 is unchanged, 2 is negative)
  for (i in 1:length(indexes)){
    
      
    # For permuting rows
     distsThisRow <- t(apply(Y, MARGIN=1, FUN=function(x) c(sum((x-X[i,])^2), sum((-1*x-X[i,])^2 ))))
     allDists[i,,] <- distsThisRow
    
  }
  
  ### To find the desired distance, we will solve the assignment problem for each pair of matrices (Y, X_{i})
  ### where X_{i} is a permuted, signed version of X and 1 <= i <= 2^{p}.
  ### The needed cost matrix can be constructed from the above distance tensor 
  ### This ensures that no row of X gets picked twice (since we don't want both signed versions of the same row being chosen)
  allmatchs <- expand.grid(rep(list(1:2), dim(X)[1]))
  
  distsListified <- list()
  
  ### Iterate over all possible variations of X
  for (row in 1:dim(allmatchs)[1]){
    
    matrixThisRow <- array(NA, dim=rep(dim(allmatchs)[2], 2))
    
    ### We've chosen the version of X to look it, so now we form the appropriate distance matrix.
    for (i in 1:length(allmatchs[row, ])){
      
      matrixThisRow[, i] <- allDists[, i, allmatchs[row, i]]
      
    }
    
    distsListified <- c(distsListified, list(matrixThisRow))
    
  } 
  
  ### Solving the allocation problem for each distance matrix.
  allPossibleSolutions <- lapply(distsListified, FUN=function(x) solve_LSAP(x))
  
  ### Calculating the distance values from the assignment rules
  potentialsolVals <- c()
  
  for (s in 1:length(allPossibleSolutions)){
  
    solValVec <- sapply(1:dim(X)[1], FUN = function(i) distsListified[[s]][i, allPossibleSolutions[[s]][i]])
    
    potentialsolVals <- c(potentialsolVals, sum(solValVec))
  
  }
  ### Return the minimal cost
  return(min(potentialsolVals))
  
}


source("Function Data Simulation Parameter Makers.R") ### Sourcing external code files. 
source("Tweedie's Formula Multivariate and Score Function.R")


literal_eval <- function(input, allowed = c("list", "c", "+", "-", "/", "*", "rep")) {
  # Create safe empty environment
  safe_env <- new.env(parent = emptyenv())
  
  # assign allowed functions
  lapply(allowed,function(f) assign(f,get(f, "package:base"),safe_env))
  
  # Evaluate input
  safe_env$expr <- parse(text = input)
  eval(substitute(expr,env = safe_env), env = safe_env)
}

twoEntriesToValParenth <- function(twoVals){ ### For formatting distance tables to LATEX
  
  paste(twoVals[1]," (",twoVals[2],")", sep="")
  
}

sim_params <- read.xlsx("functional_simulation_parameters.xlsx", sheet="Score Function Simulations",
                        rowNames = T) ### Read in simulation Parameters

## The distribution for each coordinate of the unmixed data
         
onePointFiveNorm <- function(){
    
    return(1.5)
}         
         
unmixedDist <- get(sim_params["distribution",])

### Associated parameters for the above distribution (accepted as a list)
distributionParameters <- literal_eval(sim_params['parameters',])

### Number of simulation iterations
simIterations <- as.numeric(sim_params['iterations', ])

### Sample size
n_curves <- as.numeric(sim_params['n_curves', ])
         
num_restarts <- as.numeric(sim_params['num_restarts', ])
         

### This sets the correlation parameter for mu. Set to zero for independent coordinates
rho_mu <- as.numeric(sim_params["rho_mu", ])

### Creates Sigma_mu
         
takeAbs <- as.logical(sim_params["take_abs", ])
         
makeSymmetric <- as.logical(sim_params["make_symmetric", ])

### Dimension of the vectors
p <- as.numeric(sim_params['p', ])

### Parameters for our score function estimation and ICA methods
maxIterations <- 100

gridStep = .001

lambdaGrid=10^seq(-6, 6, 1)

numberOfFolds = 7

numberOfKnotsScoreFunction = 8

nFolds <- 5

Sigma_mu <- make_rho_mat(rho_mu, p)

algorithmSpecs=c(min_iterations=as.numeric(sim_params["min_iterations_algorithm", ]), 
                 max_iterations=as.numeric(sim_params["max_iterations_algorithm", ]), 
                 min_change=as.numeric(sim_params["min_change_algorithm", ]))
  
updateUSpecs=c(min_iterations=as.numeric(sim_params["min_iterations_Uupdate", ]), 
                 max_iterations=as.numeric(sim_params["max_iterations_Uupdate", ]), 
                 min_change=as.numeric(sim_params["min_change_Uupdate", ]))   

### This is where you can experiment with different priors.
### I would suggest settng the variables "takeAbs" and "makeSymmetric"
### equal to FALSE when experimenting. 

### This code defines the quantile function for the coordinates of the unmixed prior (denoted omega).
#### To specify a prior, pass the function used for generating random samples (eg rnorm, rgamma, etc.)
#### to "do.call", along with a list (not a numeric vector, a named list) of relevant parameters, called "distribution parameters".
### For an example of what this looks like, take a look at the parameter generation file. The quantity "parameters" is an example of a named list.

### What's happening here is that we first take a large sample from the desired distribution
### and from that sample calculate the quantile function, with qarb

normSample <- rnorm(10000)         
# do.call(unmixedDist, c(list(n=10000), distributionParameters))         

   ### Good things: chi-square with df=2, gamma with shape=1, scale=5     

### This section of code samples from the desired distribution
if (takeAbs){         

    bigSample <- abs(do.call(unmixedDist, c(list(n=100000), distributionParameters)))
    
    }else if (makeSymmetric){
    
    radOnes <- sample(c(-1,1), size=100000, replace=T, prob=c(.5, .5))
    
    bigSample <- radOnes*do.call(unmixedDist, c(list(n=100000), distributionParameters))
    
    
}else{
    
    bigSample <- do.call(unmixedDist, c(list(n=100000), distributionParameters))
    
}



#clusterIndicator <- sample(c(0,1), size=100000, replace=T)

#bigSample <- rgamma(100000, shape=1, scale=.5)*clusterIndicator +(rgamma(100000, shape=1, scale=8)+10)*(1-clusterIndicator)

### Generates the quantile function from the sample
#theQuantileFunction <- qarb(bigSample) 

coordinatePic <- ggplot(NULL, aes(x=bigSample)) + geom_histogram(aes(y=..density..), bins=30, color='blue') + theme_bw(base_size=20) + 
xlab('Omega') + ylab("Density") + ggtitle("Prior Unmixed Coordinate Example")

# c(.1, .25, .5, 1, 4, 9, 16, 25, 36)  
SNRs <- seq(.5, 4, length.out=5) ### Various SNR's for which to calculate risk, L2-norm, etc.

methodNames <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}", 'SMOOTHED')   
         

### Score function assessment data
performanceDataBySNRScoreFunction <- array(NA, dim=c(length(methodNames[1:4]), simIterations, length(SNRs)))         
         
dimnames(performanceDataBySNRScoreFunction)[[1]] <- methodNames[1:4]
         
dimnames(performanceDataBySNRScoreFunction)[[2]] <- paste('itertation', 1:simIterations, sep='_')        
         
dimnames(performanceDataBySNRScoreFunction)[[3]] <- SNRs  



### Curve estimation assessment data         
performanceDataBySNRCurves <- array(NA, dim=c(length(methodNames), simIterations, length(SNRs)))         
         
dimnames(performanceDataBySNRCurves)[[1]] <- methodNames
         
dimnames(performanceDataBySNRCurves)[[2]] <- paste('itertation', 1:simIterations, sep='_')             
         
dimnames(performanceDataBySNRCurves)[[3]] <- SNRs            
         
         
### Relative risk assessment data
performanceDataBySNRRelativeRisk <- array(NA, dim=c(length(methodNames), simIterations, length(SNRs)))         
         
dimnames(performanceDataBySNRRelativeRisk)[[1]] <- methodNames
         
dimnames(performanceDataBySNRRelativeRisk)[[2]] <- paste('itertation', 1:simIterations, sep='_')             
         
dimnames(performanceDataBySNRRelativeRisk)[[3]] <- SNRs           
                

initUnmixDistFromTruth <- c()

finalUnmixDistFromTruth <- c()
         

for (snr in SNRs){    
    
    ### Initializing error arrays for a particular SNR

    errorDataScoreFunction <- array(NA, dim=c(length(methodNames[1:4]), 1, simIterations))

    errorDataCurves <- array(NA, dim=c(length(methodNames), 1, simIterations))
    
    relativeRiskDataCurves <- array(NA, dim=c(length(methodNames), 1, simIterations))
    
    dimnames(errorDataScoreFunction)[[1]] <- methodNames[1:4]

    dimnames(errorDataCurves)[[1]] <- methodNames      
    
    dimnames(relativeRiskDataCurves)[[1]] <- methodNames
    
    
    dimnames(errorDataScoreFunction)[[2]] <- list('Score Function Risk')

    dimnames(errorDataCurves)[[2]] <- list('Curve Risk')
              
    dimnames(relativeRiskDataCurves)[[2]] <- list('Relative Curve Risk')

        for (i in 1:simIterations){

          ### Iteration Start
          timeStart <- Sys.time() ### Keeping track of how long things take to run

          simulationData <- curve_generator_forceICA(n_obs=n_curves, Sigma_mu=Sigma_mu, 
                                                    SNR=snr, sigma_e =0,  ### Generating data
                                                    unmixedCoordinateDist=theQuantileFunction,
                                                   times=seq(.01, 1, .01))  

          XUnmixed <- simulationData$Wtheta %*% simulationData$theta ### Initializing the data

          X <- simulationData$theta

          trueW <- simulationData$Wtheta

          Sigma_gamma <- simulationData$Sigma_gamma

          modelingBasis <- 10*simulationData$S


          covThetaEst <- cov(t(X))

            ### If not using tweedie setting, set modelingBasis <- diag(rep(1, dim(Sigma_gamma)[1]))

          oracleTweedieInfo <- tweediesFormulaOracle(X=X, W=trueW, U=simulationData$U, Sigma_gamma=Sigma_gamma, 
                                                  functionalBasis=modelingBasis, 
                                             gridStep = .001, lambdaGrid=10^seq(-6, 3, 1), numberOfFolds = 2, 
                                                  numberOfKnotsScoreFunction = 8, maxIterations=100)
            
          oracleTweedieValues <- oracleTweedieInfo$tweedieEstimates

          ### Score function estimation assuming coordinates of X are independent

          assumedIndependentInfo <- tweedieCorrectionNonJointICA(X=X, Sigma_gamma=Sigma_gamma,  ### Tweedie without joint unmixing and score function estimation
                                             functionalBasis=modelingBasis,  
                                             transformation="none",
                                             gridStep = .001, 
                                             lambdaGrid=10^seq(-6, 3, 1),
                                             numberOfFolds = 2, 
                                             numberOfKnotsScoreFunction = 8)  


          scoreValuesAssumedIndependent <- assumedIndependentInfo$tweedieEstimates

          ### Score function estimation assuming decorrelated data is independent data

          uncorrelatedAssumedIndependentInfo <- tweedieCorrectionNonJointICA(X=X, 
                                             Sigma_gamma=Sigma_gamma,  
                                             functionalBasis=modelingBasis,  
                                             transformation="decorrelate",
                                             gridStep = .001, 
                                             lambdaGrid=10^seq(-6, 3, 1),
                                             numberOfFolds = 2, 
                                             numberOfKnotsScoreFunction = 8)

          scoreValuesUncorrAssumedIndependent <- uncorrelatedAssumedIndependentInfo$tweedieEstimates


          tweedieAlternateICAInfo <-   tweedieCorrectionNonJointICA(X=X, 
                                             Sigma_gamma=Sigma_gamma,  ### Tweedie without joint unmixing and score function estimation
                                             functionalBasis=modelingBasis,  
                                             transformation="ica", gridStep = .001, 
                                             lambdaGrid=10^seq(-6, 3, 1), numberOfFolds = 2, 
                                             numberOfKnotsScoreFunction = 8)

          scoreValuesAlternateICA <- tweedieAlternateICAInfo$tweedieEstimates

            ourMethodResults <- tweedieCorrectionWithICAWarmStart(X, Sigma_gamma=Sigma_gamma,  ### Tweedie's formula with estimtated score function.
                                           functionalBasis=modelingBasis,        ### Score function estimated with ICA approach
                                           gridStep = .001, lambdaGrid=10^seq(-6, 3, 1),
                                           numberOfFolds = 2, numberOfKnotsScoreFunction = 8,
                                        algorithmSpecs=algorithmSpecs,
                                          updateUSpecs=updateUSpecs,
                                                                 learningRateU=.2)  
            
            
  #        ourMethodResults <- tweedieCorrectionWithICARandomRestart(X, Sigma_gamma,  ### Tweedie's formula with estimtated score function.
  #                                   functionalBasis=modelingBasis,        ### Score function estimated with ICA approach
  #                                   gridStep = .001, 
  #                                   lambdaGrid=10^seq(-6, 6, 1),
  #                                   numberOfFolds = 5, 
  #                                   numberOfKnotsScoreFunction = 8,
  #                                   numRestarts=5,
  #                                   numFoldsRestarts=5,
  #                                algorithmSpecs=algorithmSpecs,
  #                                updateUSpecs=updateUSpecs)



          scoreValuesOurMethod <- ourMethodResults$tweedieEstimates

          finalW <- ourMethodResults$W

          ### Score function estimation performance
          
          integralMatrix <- (t(simulationData$S) %*% simulationData$S)/dim(simulationData$S)[1]
            
          integralMatrixOneHalf <- sqrtm(integralMatrix) ### Needed for calculating score function Risk
            
          scoreValuesThisIteration <- abind(integralMatrixOneHalf %*% scoreValuesAssumedIndependent, 
                                                      integralMatrixOneHalf %*% scoreValuesUncorrAssumedIndependent, 
                                                      integralMatrixOneHalf %*% scoreValuesAlternateICA, 
                                                      integralMatrixOneHalf %*%scoreValuesOurMethod, along=3)
            
          curvesThisIteration <- abind(t(simulationData$S %*% scoreValuesAssumedIndependent), t(simulationData$S %*%
                                            scoreValuesUncorrAssumedIndependent), 
                                            t(simulationData$S %*% scoreValuesAlternateICA), 
                                            t(simulationData$S %*% scoreValuesOurMethod),
                                       t(simulationData$S %*% X), along=3 ) 
            
            
          oracleCurvesThisIteration <- t(simulationData$S %*% oracleTweedieValues )
            
            
            
          #### Score Function Risk
          errorsThisIterationScoreFunction <- plyr::aaply(scoreValuesThisIteration, .margins=3, 
                                          .fun= averagedVectorl2Norm, Y=integralMatrixOneHalf%*% oracleTweedieValues)
            
          aggErrorsThisIterationScoreFunction <- rowMeans(errorsThisIterationScoreFunction)  
            
          ### Error, Tweedie estimators to truth
          errorsThisIterationCurves <- plyr::aaply(curvesThisIteration, .margins=3, 
                                          .fun= averagedVectorl2Norm, 
                                                      Y=t(simulationData$S %*% simulationData$mu))
          
          aggErrorsThisIterationCurves <- rowMeans(errorsThisIterationCurves)   
            
            
            ### Error, oracle tweedie to truth  
          riskTweedieToTruth <- errorsThisIterationCurves
            
          riskOracleToTruth <- averagedVectorl2Norm(X=oracleCurvesThisIteration, Y=t(simulationData$S %*% simulationData$mu))
            
          riskOracleToTruth <- sapply(riskOracleToTruth, FUN=rep, each=length(methodNames))
            
          aggRelativeRisksThisIteration <- rowMeans((riskTweedieToTruth-riskOracleToTruth)/riskOracleToTruth)
            
            
            
          ### Populating error data arrays with the appropriate values  
          performanceDataBySNRScoreFunction[, i, as.character(snr)] <- aggErrorsThisIterationScoreFunction

          performanceDataBySNRCurves[, i, as.character(snr)] <- aggErrorsThisIterationCurves
            
          performanceDataBySNRRelativeRisk[, i, as.character(snr)] <-  aggRelativeRisksThisIteration 

          ### Unmixing matrix estimation performance
          initUnmixDistFromTruth <- c(initUnmixDistFromTruth, distWithPermutAndSgnChange(X=tweedieAlternateICAInfo$W, Y=trueW))

          finalUnmixDistFromTruth <- c(finalUnmixDistFromTruth, distWithPermutAndSgnChange(X=finalW, Y=trueW ))

          timeEnd <- Sys.time()

          outMessage <- paste("SNR: ", snr, '; Iteration: ', i, sep='')
      
          print(outMessage)
            
          print(difftime(timeEnd, timeStart, units='mins'))

        }
                                        
}
         
         

 ### Score func MSE against SNR        
meltedPerformanceDataBySNRScoreFunction <- melt(performanceDataBySNRScoreFunction)
names(meltedPerformanceDataBySNRScoreFunction) <- c('Method', 'Iteration', 'SNR', 'Risk')

meltedPerformanceDataBySNRScoreFunction <- meltedPerformanceDataBySNRScoreFunction %>% filter(!is.na(Risk))

aggMeltedPerformanceDataBySNRScoreFunction <- meltedPerformanceDataBySNRScoreFunction %>% group_by(Method, SNR=round(SNR, 3)) %>% 
summarize(SE=sciFormatter(sd(Risk, na.rm=T)/sqrt(n()), 3)   , Risk=sciFormatter(max(0, mean(Risk, na.rm=T)), 3),
         `Risk (SE)` = paste(Risk, ' (', SE, ')', sep=''))
         
meltedPerformanceDataBySNRScoreFunction$Method <-factor(
    meltedPerformanceDataBySNRScoreFunction$Method, levels=methodNames)


snrScoreFuncPlot <- ggplot(aggMeltedPerformanceDataBySNRScoreFunction, aes(x=as.numeric(SNR), y=as.numeric(Risk), color=Method)) + geom_line(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('MSE') + ggtitle("MSE for Score Function Against SNR, by Method")

snrScoreFuncPlotWithUncertainty <- ggplot(meltedPerformanceDataBySNRScoreFunction, aes(x=SNR, y=Risk, color=Method)) + geom_smooth(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('MSE') + ggtitle("MSE for Score Function Against SNR, by Method")



### Score function error table, formatted for latex
aggMeltedPerformanceDataBySNRScoreFunction <- reshape2::dcast(aggMeltedPerformanceDataBySNRScoreFunction, Method ~ SNR, value.var="Risk (SE)")

aggMeltedPerformanceDataBySNRScoreFunctionForLatex <- knitr::kable(aggMeltedPerformanceDataBySNRScoreFunction, 'latex')

### Relative Risk against SNR        
meltedPerformanceDataBySNRRelativeRisk<- melt(performanceDataBySNRRelativeRisk)
         
names(meltedPerformanceDataBySNRRelativeRisk) <- c('Method', 'Iteration', 'SNR', 'Relative Risk')

meltedPerformanceDataBySNRRelativeRisk <- meltedPerformanceDataBySNRRelativeRisk %>% filter(!is.na(`Relative Risk`))
         
meltedPerformanceDataBySNRRelativeRisk$Method <-factor(
    meltedPerformanceDataBySNRRelativeRisk$Method, levels=methodNames)


aggMeltedPerformanceDataBySNRRelativeRisk <- meltedPerformanceDataBySNRRelativeRisk%>% group_by(Method, SNR=round(SNR, 3)) %>% 
summarize(SE=sciFormatter(sd(`Relative Risk`, na.rm=T)/sqrt(n()), 3), `Relative Risk`=sciFormatter( mean(`Relative Risk`, na.rm=T), 3),
         `Relative Risk (SE)` = paste(`Relative Risk`, ' (', SE, ')', sep=''))



snrRelativeRiskPlot <- ggplot(aggMeltedPerformanceDataBySNRRelativeRisk, aes(x=as.numeric(SNR), y=as.numeric(`Relative Risk`), color=Method)) +
geom_line(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('MSE') + ggtitle("Relative Risk Against SNR, by Method")

snrRelativeRiskPlotWithUncertainty <- ggplot(meltedPerformanceDataBySNRRelativeRisk, aes(x=SNR, y=( `Relative Risk`), color=Method)) +
geom_smooth(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('MSE') + ggtitle("Relative Risk Against SNR, by Method")
         

aggMeltedPerformanceDataBySNRRelativeRisk <- reshape2::dcast(aggMeltedPerformanceDataBySNRRelativeRisk, Method ~ SNR, value.var="Relative Risk (SE)")

aggMeltedPerformanceDataBySNRRelativeRiskForLatex <- knitr::kable(aggMeltedPerformanceDataBySNRRelativeRisk, 'latex')

#### Creates the file directory for storing the results. 


rootOutputDir <- "../reports/score_and_curve_snr_reports"                                                                              
                                                                              
if(!dir.exists(rootOutputDir)){
  
  dir.create(rootOutputDir, recursive=TRUE)
  
  
}      
         
         
distributionDirectory <- str_trim(paste(c("", "Abs")[(takeAbs ==TRUE)+1], paste(c("", "Symmetrized")[(makeSymmetric ==TRUE)+1],
    c("", "Mixture")[any(lengths(distributionParameters) > 1)+1], 
                               str_replace(sim_params["distribution",], "^r", ""))))


nodeName <- Sys.info()["nodename"]

rootToDistribution <- paste(rootOutputDir, nodeName, distributionDirectory, sep='/')

if(!dir.exists(rootToDistribution)){
  
  dir.create(rootToDistribution, recursive = T)
  
  
}

runDate <- Sys.Date()

directoryFileNames <- list.files(rootToDistribution)

howManyThisDateSoFar <- sum(str_detect(directoryFileNames, pattern=as.character(runDate)))

outputDirectory <- paste(rootToDistribution, paste(Sys.Date(), ": ",howManyThisDateSoFar+1, sep=''), sep='/')

dir.create(outputDirectory)       

### Results from different analyses go to different directories. 
scoreDirectory <- paste(outputDirectory, 'score_function_estimation_performance', sep='/')                                                                              

relativeRiskDirectory <- paste(outputDirectory, 'relative_risk_performance', sep='/')  
                                                                              
if(!dir.exists(scoreDirectory)){
  
  dir.create(scoreDirectory)
  
  
}                
        
if(!dir.exists(relativeRiskDirectory )){
  
  dir.create(relativeRiskDirectory )
  
  
}                           
  



#### Outputting Score function info

                                                                              
ggsave(plot=snrScoreFuncPlot, filename=paste(scoreDirectory, 
        "score_function_estimation_against_snr.pdf", sep='/'), 
       units='in', width=8, height=8)
         
ggsave(plot=snrScoreFuncPlotWithUncertainty, filename=paste(scoreDirectory, 
        "score_function_estimation_against_snr_uncertainty.pdf", sep='/'), 
       units='in', width=8, height=8)    
         
fileConn<-file(paste(scoreDirectory, 
        "score_against_snr_table.txt", sep='/'))
writeLines(tableFixer(aggMeltedPerformanceDataBySNRScoreFunctionForLatex)  , fileConn)
close(fileConn)       

### Outputting relative risk estimation info


ggsave(plot=snrRelativeRiskPlot, filename=paste(relativeRiskDirectory, 
        "relative_risk_against_snr.pdf", sep='/'), 
       units='in', width=8, height=8) 
         
         
ggsave(plot=snrRelativeRiskPlotWithUncertainty, filename=paste(relativeRiskDirectory, 
        "relative_risk_against_snr_uncertainty.pdf", sep='/'), 
       units='in', width=8, height=8)

fileConn<-file(paste(relativeRiskDirectory, 
        "relative_risk_against_snr_table.txt", sep='/'))
writeLines(tableFixer(aggMeltedPerformanceDataBySNRRelativeRiskForLatex)  , fileConn)
close(fileConn)     

### Outputting simulation parameters.

write.csv(sim_params, file=paste(outputDirectory, 'parameters.csv', sep='/'))



### Example of unmized coordinate distribution, sans Gaussian noise
ggsave(plot=coordinatePic, filename=paste(outputDirectory, 
        "unmixed_coordinate_distribution.pdf", sep='/'), 
       units='in', width=8, height=8)

stopCluster(cl)


