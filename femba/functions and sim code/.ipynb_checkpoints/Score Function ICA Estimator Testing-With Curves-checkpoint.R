### Good Sims: gamma(shape=1, scale=5)
library(stringr)
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

averagedVectorl2Norm <- function(X,Y){

    return(rowMeans((X-Y)^2))
    
}

source("Function Data Simulation Parameter Makers.R")
source("Tweedie's Formula Multivariate and Score Function.R")


desired_clusters <- 8

max_available_clusters <- detectCores()

cl <- makeCluster(min(c(max_available_clusters, desired_clusters))-3)

registerDoParallel(cl)

scoreListToVector <- function(aScoreList){
  
  scoreVector <- function(theta, derivativeOrder=0){
    
    if (!is.matrix(theta)){
      
      theta <- array(theta, dim=c(length(theta), 1))
      
    }
    
    t(sapply(1:dim(theta)[1], FUN=function(x) aScoreList[[x]](theta[x,], derivativeOrder=derivativeOrder)))
    
  }
  
  return(scoreVector)
  
  
}

### Calculates minimum distance between a matrix Y, and the row-permuted, signed version of a matrix, X. 
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

### Data generator for ICA simulations (quantile_func is a quantile function)
ICAmodelDataMaker <- function(n, p, quantile_func, rho=.5, rho_misspecification=0, rho_gamma=.2){

  misspecificationMat <- make_rho_mat(rho_misspecification, p)
  
  #XUnmixed <- matrixToPower(misspecificationMat, -.5) %*% XUnmixedPure
  
  ### True U in so(p)
  UTrue <- gramSchmidt(array(rnorm(p^2), dim=rep(p, 2)))$Q 
  
  ### True covariance matrix for data
  covMat <- make_rho_mat(rho, dim(UTrue)[1])
  
  Sigma_g <- make_rho_mat(rho_gamma, dim(UTrue)[1])
  
  ### Always the identity matrix for these simulations
  sigmaTilde <- Sigma_g %*% Sigma_g
  
  ### Covariance matrix to the -1/2 power
  covMatMinusHalf <- matrixToPower(covMat, -.5)
  
  ### Sigma tilde to the -1/2 power
  sigmaTildeMinusHalf <- matrixToPower(sigmaTilde, -.5)
    
  ### True unmixing matrix
  trueW <- UTrue %*% covMatMinusHalf
  
  ### True unmixing matrix for data with covariance Sigma tilde
  trueWTilde <- UTrue %*% sigmaTildeMinusHalf
    
  ### Generate unmixed data with NORTA
  XUnmixed <- norta(n, corr_mat = misspecificationMat, distribution = quantile_func)
  
  XUnmixed <- matrixToPower(cov(t(XUnmixed)), -.5)  %*% XUnmixed
    
  XUnmixed <- XUnmixed - t(sapply(rowMeans(XUnmixed), rep, each=dim(XUnmixed)[2]))
    
  ### The observed, mixed data
  X <- solve(trueW) %*% XUnmixed 
  
  oneVec <- matrix(1, nrow=dim(UTrue)[2], ncol=1)
  
  ### The true constants needed for estimating the score function
  cTrueWTilde <- apply(trueWTilde , MARGIN=1, 
                       FUN = function(x) ( t(x) %*% sigmaTilde %*% t(t(x))))
  
  finalOutput <- list(XUnmixed=XUnmixed, trueW=trueW, trueWTilde=trueWTilde, 
                      sigmaTilde=sigmaTilde, X=X, cTrueWTilde=cTrueWTilde)

}

### Similar to ast.literal eval in Python
### Used for data simulation from input datailing what distribution data should come from.
### See parameter file for more details.
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

### Read in parameter file
sim_params <- read.xlsx("functional_simulation_parameters.xlsx", sheet="Score Function Simulations",
                        rowNames = T)

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
         
rho_misspecification <- as.numeric(sim_params["rho_misspecification", ])
         
         
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
         
         
#STOPPING CRITERIA FOR ALGORITHM
algorithmSpecs=c(min_iterations=as.numeric(sim_params["min_iterations_algorithm", ]), 
                 max_iterations=as.numeric(sim_params["max_iterations_algorithm", ]), 
                 min_change=as.numeric(sim_params["min_change_algorithm", ]))
  
updateUSpecs=c(min_iterations=as.numeric(sim_params["min_iterations_Uupdate", ]), 
                 max_iterations=as.numeric(sim_params["max_iterations_Uupdate", ]), 
                 min_change=as.numeric(sim_params["min_change_Uupdate", ]))   

normSample <- rnorm(10000)         
# do.call(unmixedDist, c(list(n=10000), distributionParameters))         

   ### Good things: chi-square with df=2, gamma with shape=1, scale=5     
if (takeAbs){         

    newThing <- abs(do.call(unmixedDist, c(list(n=100000), distributionParameters)))
    
    }else if (makeSymmetric){
    
    radOnes <- sample(c(-1,1), size=100000, replace=T, prob=c(.5, .5))
    
    newThing <- radOnes*do.call(unmixedDist, c(list(n=100000), distributionParameters))
    
    
}else{
    
    newThing <- do.call(unmixedDist, c(list(n=100000), distributionParameters))
    
}
## Try a skew normal distribution, where f(x)=2phi(x)Phi(alpha*x)  
## So, 
         
#newThing <- abs(rnorm(100000))^1.5
         
theQuantileFunction <- qarb(newThing)

initUnmixDistFromTruth <- c()

finalUnmixDistFromTruth <- c()
         

         
         
SNRs <- c(.1, .25, .5, 1, 4, 9, 16, 25, 36)         
       
         
performanceDataBySNRScoreFunction <- array(NA, dim=c(4, 1, length(SNRs)))         
         
dimnames(performanceDataBySNRScoreFunction)[[1]] <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}") 
         
dimnames(performanceDataBySNRScoreFunction)[[2]] <- c("risk")          
         
dimnames(performanceDataBySNRScoreFunction)[[3]] <- SNRs   
         
performanceDataBySNRCurves <- array(NA, dim=c(4, 1, length(SNRs)))         
         
dimnames(performanceDataBySNRCurves)[[1]] <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}") 
         
dimnames(performanceDataBySNRCurves)[[2]] <- c("l2_norm")          
         
dimnames(performanceDataBySNRCurves)[[3]] <- SNRs            
         
         
performanceDataBySNRRelativeRisk <- array(NA, dim=c(4, 1, length(SNRs)))         
         
dimnames(performanceDataBySNRRelativeRisk)[[1]] <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}") 
         
dimnames(performanceDataBySNRRelativeRisk)[[2]] <- c("l2_norm")          
         
dimnames(performanceDataBySNRRelativeRisk)[[3]] <- SNRs           
         
         
         
methodNames <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}")          
         
for (snr in SNRs){         

    errorDataScoreFunction <- array(NA, dim=c(4, 1, simIterations))

    errorDataCurves <- array(NA, dim=c(4, 1, simIterations))
    
    relativeRiskDataCurves <- array(NA, dim=c(4, 1, simIterations))
    
    
    

    dimnames(errorDataScoreFunction)[[1]] <- methodNames

    dimnames(errorDataCurves)[[1]] <- methodNames      
    
    dimnames(relativeRiskDataCurves)[[1]] <- methodNames
    
    
    dimnames(errorDataScoreFunction)[[2]] <- list('Score Function Risk')

    dimnames(errorDataCurves)[[2]] <- list('Curve Risk')
              
    dimnames(relativeRiskDataCurves)[[2]] <- list('Relative Curve Risk')

        for (i in 1:simIterations){

          ### Iteration Start
          timeStart <- Sys.time()

          simulationData <- curve_generator_forceICA(n_obs=n_curves, Sigma_mu=make_rho_mat(rho=.5, p),
                                                    SNR=snr, sigma_e =0, 
                                                    unmixedCoordinateDist=theQuantileFunction,
                                                   times=seq(.01, 1, .01))  

          XUnmixed <- simulationData$Wtheta %*% simulationData$theta

          X <- simulationData$theta

          trueW <- simulationData$Wtheta

          Sigma_gamma <- simulationData$Sigma_gamma

          modelingBasis <- 10*simulationData$S


          covThetaEst <- cov(t(X))

            ### If not using tweedie setting, set modelingBasis <- diag(rep(1, dim(Sigma_gamma)[1]))

          covThetaMinusHalfEst <- matrixToPower(covThetaEst, -.5)

          trueTweedieCorrectionInfo <- tweediesFormulaOracle(X=X, W=trueW, U=trueW %*% solve(covThetaMinusHalfEst), Sigma_gamma=Sigma_gamma, 
                                                          functionalBasis=modelingBasis, 
                                                     gridStep = .001, lambdaGrid=10^seq(-6, 3, 1), numberOfFolds = 2, 
                                                          numberOfKnotsScoreFunction = 8, maxIterations=100)

          oracleTweedieValues <- trueTweedieCorrectionInfo$tweedieEstimates

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

          ### Score function estimation with our method
         # ourMethodResults <- tweedieCorrectionWithICARandomRestarts(X, Sigma_gamma=Sigma_gamma,  ### Tweedie's formula with estimtated score function.
             #                                functionalBasis=diag(rep(1, dim(Sigma_gamma)[1])),        ### Score function estimated with ICA approach
              #                               gridStep = .001, 
            ##                                 lambdaGrid=10^seq(-6, 6, 1),
            #                                 numberOfFolds = 5, 
            #                                 numberOfKnotsScoreFunction = 8,
            #                                 numRestarts=num_restarts,
            ##                                 numFoldsRestarts=2,
            #                              algorithmSpecs=c(min_iterations=1, max_iterations=100, min_change=1e-06),
            #                              updateUSpecs=c(min_iterations=1, max_iterations=1, min_change=1e-06))




            ourMethodResults <- tweedieCorrectionWithICAWarmStart(X, Sigma_gamma=Sigma_gamma,  ### Tweedie's formula with estimtated score function.
                                             functionalBasis=modelingBasis,        ### Score function estimated with ICA approach
                                             gridStep = .001, lambdaGrid=10^seq(-6, 3, 1),
                                             numberOfFolds = 2, numberOfKnotsScoreFunction = 8,
                                          algorithmSpecs=algorithmSpecs,
                                          updateUSpecs=updateUSpecs,
                                                                 learningRateU=.2)  



          scoreValuesOurMethod <- ourMethodResults$tweedieEstimates

          finalW <- ourMethodResults$W

          ### Score function estimation performance
          
          integralMatrix <- (t(simulationData$S) %*% simulationData$S)/dim(simulationData$S)[1]
            
          integralMatrixOneHalf <-  matrixToPower(integralMatrix, .5)
            
          scoreValuesThisIteration <- abind(integralMatrixOneHalf %*% scoreValuesAssumedIndependent, 
                                                      integralMatrixOneHalf %*% scoreValuesUncorrAssumedIndependent, 
                                                      integralMatrixOneHalf %*% scoreValuesAlternateICA, 
                                                      integralMatrixOneHalf %*%scoreValuesOurMethod, along=3)
            
          curvesThisIteration <- abind(t(simulationData$S %*% scoreValuesAssumedIndependent), t(simulationData$S %*%
                                            scoreValuesUncorrAssumedIndependent), 
                                            t(simulationData$S %*% scoreValuesAlternateICA), 
                                            t(simulationData$S %*% scoreValuesOurMethod), along=3 ) 
            
            
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
          riskTweedieToOracle <- errorsThisIterationCurves
            
          riskOracleToTruth <- averagedVectorl2Norm(X=oracleCurvesThisIteration, Y=t(simulationData$S %*% simulationData$mu))
            
          riskOracleToTruth <- sapply(riskOracleToTruth, FUN=rep, each=4)
            
          aggRelativeRisksThisIteration <- rowMeans((riskTweedieToOracle-riskOracleToTruth)/riskOracleToTruth)  

          errorDataScoreFunction[, , i] <- aggErrorsThisIterationScoreFunction

          errorDataCurves[, , i] <- aggErrorsThisIterationCurves
            
          relativeRiskDataCurves[, , i] <- aggRelativeRisksThisIteration 

          ### Unmixing matrix estimation performance
          initUnmixDistFromTruth <- c(initUnmixDistFromTruth, distWithPermutAndSgnChange(X=tweedieAlternateICAInfo$W, Y=trueW))

          finalUnmixDistFromTruth <- c(finalUnmixDistFromTruth, distWithPermutAndSgnChange(X=finalW, Y=trueW ))

          timeEnd <- Sys.time()

          print(paste("Iteration", i))
          print(difftime(timeEnd, timeStart, units='mins'))

        }
    
    
    
    
    
    
    

    ### Aggregating simulation error data      
    
    
    
    #### Score Function Risk
    averageTableScoreFunction <- plyr::aaply(errorDataScoreFunction, .margins=c(1,2), .fun=mean)

    performanceDataScoreFunction <- data.frame(errorDataScoreFunction)

    performanceDataScoreFunction$Method <- row.names(performanceDataScoreFunction)

    performanceDataScoreFunction <- melt(setDT(performanceDataScoreFunction), id.vars = "Method", variable.name = "iteration",
                           value.name="Risk_that_iteration")[,c(1,3)]         


    ### Curve Risk
    averageTableCurves <- plyr::aaply(errorDataCurves, .margins=c(1,2), .fun=mean)

    performanceDataCurves <- data.frame(errorDataCurves)

    performanceDataCurves$Method <- row.names(performanceDataCurves)

    performanceDataCurves <- melt(setDT(performanceDataCurves), id.vars = "Method", variable.name = "iteration",
                           value.name="L2_norm_that_iteration")[,c(1,3)]             


    ### Relative Curve Risk
    averageTableRelativeCurveRisk <- plyr::aaply(relativeRiskDataCurves, .margins=c(1,2), .fun=mean)

    performanceDataRelativeRisk<- data.frame(relativeRiskDataCurves)

    performanceDataRelativeRisk$Method <- row.names(relativeRiskDataCurves)

    performanceDataRelativeRisk <- melt(setDT(performanceDataRelativeRisk), id.vars = "Method", 
                                        variable.name = "iteration",
                           value.name="L2_norm_that_iteration")[,c(1,3)]    




    performanceDataBySNRScoreFunction[, , as.character(snr)] <-  averageTableScoreFunction                                             


    performanceDataBySNRCurves[, , as.character(snr)] <-  averageTableCurves 
    
    performanceDataBySNRRelativeRisk[, , as.character(snr)] <-  averageTableRelativeCurveRisk
    print(paste("SNR:", snr))
                                        
}
         
         
         
         
         
 ### Score func MSE against SNR        
meltedPerformanceDataBySNRScoreFunction <- melt(performanceDataBySNRScoreFunction)[, -2]
         
names(meltedPerformanceDataBySNRScoreFunction) <- c('Method', 'SNR', 'Risk')
         
meltedPerformanceDataBySNRScoreFunction$Method <-as.character(
    meltedPerformanceDataBySNRScoreFunction$Method)

snrScoreFuncPlot <- ggplot(meltedPerformanceDataBySNRScoreFunction, aes(x=SNR, y=Risk, color=Method)) + geom_line(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('MSE') + ggtitle("MSE for Score Function Against SNR, by Method")
         
         
performanceBoxplotScoreFunction <- (ggplot(performanceDataScoreFunction, aes(x=Method, y=Risk_that_iteration, fill=Method)) + 
                      geom_boxplot() + theme_bw(base_size=20) + 
                      theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()) +ylab("MSE")
                    + ggtitle("MSE for Score Function Estimation"))         
         
         
         
         
         
         
  ### Curve estimation MSE against SNR        
meltedPerformanceDataBySNRCurves <- melt(performanceDataBySNRCurves)[, -2]
         
names(meltedPerformanceDataBySNRCurves) <- c('Method', 'SNR', 'l2-norm')
         
meltedPerformanceDataBySNRCurves$Method <-as.character(
    meltedPerformanceDataBySNRCurves$Method)


### Why MSE increases as a function of SNR
### MSE between original and true data is a constant
### theta+Sigma*nu-mu= (theta-mu)+Sigma*nu
### When SNR is large, Sigma*mu is small 
snrCurvePlot <- ggplot(meltedPerformanceDataBySNRCurves, aes(x=SNR, y=sqrt(`l2-norm`), color=Method)) + 
         geom_line(lwd=1.5)  +
         theme_bw(base_size=20) + ylab('L2-Norm') + ggtitle("L2-Norm for Curves Against SNR, by Method") 
                 
         
performanceBoxplotCurves <- (ggplot(performanceDataCurves, aes(x=Method, y=sqrt(L2_norm_that_iteration), fill=Method)) + 
                      geom_boxplot() + theme_bw(base_size=20) + 
                      theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()) +ylab("L2-Norm")
                    + ggtitle("L2-Norm for Curve Estimation"))                 
         
                                                    
                                                                              
                                                                              
      
                                                                              
                                                                              
                                                                              
                                                                              
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
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
                                                                              
scoreDirectory <- paste(outputDirectory, 'score_function_estimation_performance', sep='/')                                                                              
                                                                              
curveDirectory <- paste(outputDirectory, 'curve_estimation_performance', sep='/')                                                                                
                                                                              
if(!dir.exists(scoreDirectory)){
  
  dir.create(scoreDirectory)
  
  
}
                                                                              
                                                                              
if(!dir.exists(curveDirectory )){
  
  dir.create(curveDirectory)
  
  
}                           
                                                                                                                                     
                                                                              
ggsave(plot=snrScoreFuncPlot, filename=paste(scoreDirectory, 
        "score_function_estimation_against_snr.pdf", sep='/'), 
       units='in', width=8, height=8)
         
ggsave(plot=performanceBoxplotScoreFunction, filename=paste(scoreDirectory, 
        "boxplot_performance_score_function_highest_SNR.pdf", sep='/'), 
       units='in', width=8, height=8)         
         
         
         
         
         
         
         
         
ggsave(plot=snrCurvePlot, filename=paste(curveDirectory, 
        "curve_estimation_against_snr.pdf", sep='/'), 
       units='in', width=8, height=8) 
         
         
ggsave(plot=performanceBoxplotCurves, filename=paste(curveDirectory, 
        "boxplot_performance_curves_highest_SNR.pdf", sep='/'), 
       units='in', width=8, height=8)

write.csv(sim_params, file=paste(outputDirectory, 'parameters.csv', sep='/'))

