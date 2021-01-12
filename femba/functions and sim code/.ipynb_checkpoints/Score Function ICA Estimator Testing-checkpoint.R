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

source("Function Data Simulation Parameter Makers.R")
source("Tweedie's Formula Multivariate and Score Function.R")

averagedVectorl2Norm <- function(X,Y){
    
    return((norm(X-Y,'F')^2)/dim(Y)[2])
    
}


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
         
errorData <- array(NA, dim=c(4, 6, simIterations))
    
dimnames(errorData)[[1]] <- c("\\fembant{}", "\\fembat{}", "\\fembafastICA{}", "\\fembajointICA{}") 
         
for (i in 1:simIterations){
  
  ### Iteration Start
  timeStart <- Sys.time()
  
  ### Simulate data
  simulationData <- ICAmodelDataMaker(n=n_curves, p=p, quantile_func=theQuantileFunction, rho=.5, 
                                      rho_misspecification=rho_misspecification,
                                     rho_gamma=.3)
  
  ### Extract simulation parmeters
  XUnmixed <- simulationData$XUnmixed
  
 X <- simulationData$X
  
  trueW <- simulationData$trueW
    
  Sigma_gamma <- matrixToPower(simulationData$sigmaTilde, .5)
    
  modelingBasis <- diag(rep(1, dim(Sigma_gamma)[1]))
  
  
  covThetaEst <- cov(t(X))
    
    ### If not using tweedie setting, set modelingBasis <- diag(rep(1, dim(Sigma_gamma)[1]))
  
  covThetaMinusHalfEst <- matrixToPower(covThetaEst, -.5)
  
  trueTweedieCorrectionInfo <- tweediesFormulaOracle(X=X, W=trueW, U=trueW %*% solve(covThetaMinusHalfEst), Sigma_gamma=Sigma_gamma, 
                                                  functionalBasis=modelingBasis, 
                                             gridStep = .001, lambdaGrid=10^seq(-6, 3, 1), numberOfFolds = 2, 
                                                  numberOfKnotsScoreFunction = 8, maxIterations=100)
  
  trueScoreValues <- trueTweedieCorrectionInfo$tweedieEstimates
  
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
    
    norm(scoreValuesOurMethod-trueScoreValues, "F")
    
                                                  
  finalW <- ourMethodResults$W
  
  ### Score function estimation performance
    
  scoreValuesThisIteration <- abind(scoreValuesAssumedIndependent, scoreValuesUncorrAssumedIndependent, scoreValuesAlternateICA, 
                                    scoreValuesOurMethod, along=3) 
    
  ErrorInfoThisIteration <- plyr::aaply(scoreValuesThisIteration, .margins=3, 
                                  .fun= bias_variance_mse, y=trueScoreValues)
    
  if (i==1){
      
      
      dimnames(errorData)[[2]] <- dimnames(ErrorInfoThisIteration )[[2]]
      
  }
    
  errorData[, , i] <- ErrorInfoThisIteration
  
  ### Unmixing matrix estimation performance
  initUnmixDistFromTruth <- c(initUnmixDistFromTruth, distWithPermutAndSgnChange(X=tweedieAlternateICAInfo$W, Y=trueW))
  
  finalUnmixDistFromTruth <- c(finalUnmixDistFromTruth, distWithPermutAndSgnChange(X=finalW, Y=trueW ))
    
  timeEnd <- Sys.time()
  
  print(paste("Iteration", i))
  print(difftime(timeEnd, timeStart, units='mins'))

}

         
#### DIAGNOSTICS FOR MODEL FITS 
         
## Distribution of errors for a particular coordinate
coordinateForPlotting <- 1   
         
alternateErrors <- scoreValuesAlternateICA[coordinateForPlotting ,]-trueScoreValues[coordinateForPlotting ,]
ICAErrors <- scoreValuesOurMethod[coordinateForPlotting ,]-trueScoreValues[coordinateForPlotting ,] 

ICAErrorDataNames <- c(rep("FastICA", n_curves), rep("Our ICA", n_curves))
         
ICAErrorDistrubutionData <- cbind.data.frame(ICAErrorDataNames, c(alternateErrors, ICAErrors))
         
names(ICAErrorDistrubutionData) <- c("Method", "Error")
         
coordinateErrorExample <- (ggplot(ICAErrorDistrubutionData, aes(x=Error, fill=Method))
 + geom_histogram(alpha=.5, position="identity")  
 + theme_bw(base_size=20)  + ylab("Count") + 
        ggtitle("The Error Distribution for one of the Mixed Score Function Coordinates"))
         

#### DIAGNOSTICS FOR UNMIXED DATA
WAlternateICATheta <- tweedieAlternateICAInfo$W %*% X     

unmixedAlternateICACurves <- solve(t(tweedieAlternateICAInfo$W)) %*% solve(Sigma_gamma) %*%  (scoreValuesAlternateICA-X)
         
WOurMethodTheta <- finalW %*% X   
unmixedOurICACurves <- solve(t(finalW)) %*% solve(Sigma_gamma) %*%  (scoreValuesOurMethod-X)
         
WTrueTheta <- trueW %*% X 
         
unmixedTrueCurves <- solve(t(trueW)) %*% solve(Sigma_gamma) %*%  (trueScoreValues-X)         
         
unmixedAlternateICA <- WAlternateICATheta[coordinateForPlotting, ]         
         
unmixedOurICA <- WAlternateICATheta[coordinateForPlotting, ]  

unmixedTruth <- WTrueTheta[coordinateForPlotting, ]  
         
ICAUnmixedCoordinateData <- cbind.data.frame(ICAErrorDataNames, c(unmixedAlternateICA, unmixedOurICA))
         
names(ICAUnmixedCoordinateData) <- c("Method", "Value")         
         
coordinateDistributionExample <- (ggplot(ICAUnmixedCoordinateData, aes(x=Value, fill=Method))
 + geom_histogram(alpha=.5, position="identity")  
 + theme_bw(base_size=20)   + ylab("Count") + 
    ggtitle("An Example of the Distribution for one of the Unmixed Coordinates"))         
      
         
         
         
         
         
         
         
         
         
         
         
         
         
### Comparing score functions of unmixed coordinates for each method
         
unimxedScoreFunctionCoordinateExample <- cbind.data.frame(
    c(rep("FastICA", n_curves), rep("Our ICA", n_curves), rep("Oracle", n_curves)),
    c(WAlternateICATheta[coordinateForPlotting, ],WOurMethodTheta[coordinateForPlotting, ], WTrueTheta[coordinateForPlotting, ]),
    c(unmixedAlternateICACurves[coordinateForPlotting, ],unmixedOurICACurves[coordinateForPlotting, ],
     unmixedTrueCurves[coordinateForPlotting,]) ) 
         
        
         
names(         
unimxedScoreFunctionCoordinateExample) <- c("Method", "WTheta", "Value")
         
unimxedScoreFunctionCoordinateExample <- (unimxedScoreFunctionCoordinateExample %>% arrange(Method, WTheta) )
         
onlyTruth <- unimxedScoreFunctionCoordinateExample %>% filter(Method=="Oracle")
         
onlyEstimates <- unimxedScoreFunctionCoordinateExample %>% filter(Method!="Oracle")    

onlyEstimates %>% group_by(Method) %>% summarize(MSE=mean((Value-onlyTruth$Value)^2)  )       
    
unmixedCoordinateExamplePlot <- (ggplot(unimxedScoreFunctionCoordinateExample, 
                        aes(x=WTheta, y=Value, color=Method)) + geom_smooth(method="loess")
                                +theme_bw(base_size=20) + 
        ggtitle("An Example of the Score Function for one of the Unmixed Coordinates"))     
         

         
         
         
         
         
mean((unmixedAlternateICACurves[coordinateForPlotting, ]-unmixedTrueCurves[coordinateForPlotting,])^2)         
         
         
mean((unmixedOurICACurves[coordinateForPlotting, ]-unmixedTrueCurves[coordinateForPlotting,])^2)                
         

### Aggregating simulation error data         
averageTable <- plyr::aaply(errorData, .margins=c(1,2), .fun=mean)
         
performanceData <- data.frame(errorData[, "mse", ])
         
performanceData$Method <- row.names(performanceData)
         
performanceData <- melt(setDT(performanceData), id.vars = "Method", variable.name = "iteration",
                       value.name="MSE_that_iteration")[,c(1,3)]         

         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
### Creating data frame of unmixing estimation performance values
unmixingDistValues <- c(initUnmixDistFromTruth, finalUnmixDistFromTruth)

distType <- rep(c("Standard ICA", "Our Approach"), each=length(unmixingDistValues)/2)

distancesAndTypes <- cbind.data.frame(distType, unmixingDistValues)

names(distancesAndTypes) <- c("Type", "Distance")

### Boxplot of distances
distBoxPlots <- (ggplot(distancesAndTypes, aes(x=Type, y=Distance, fill=Type)) + geom_boxplot() + theme_bw(base_size=20) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +ylab("Distance"))

### Table of distances
distanceTable <- distancesAndTypes %>% group_by(Type) %>% dplyr::summarize(`Mean Distance`=round(mean(Distance), 3), 
                                                                    SE=round(sd(Distance)/sqrt(length(simIterations)), 3))

distanceDisplayName <- twoEntriesToValParenth(names(distanceTable)[2:3])


distanceTable[, distanceDisplayName] <- apply(distanceTable[,2:dim(distanceTable)[2]], MARGIN=1, 
                                      FUN=function(x) twoEntriesToValParenth(x))

distanceForDisplay <- distanceTable[, c(1,4)]

distanceTableLatex <- kable(distanceForDisplay, "latex")

                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
                                              
### Creating table for each method's performance. 
performanceTable <- performanceData %>% group_by(Method) %>% dplyr::summarise(MSE=round(mean(MSE_that_iteration), 3), 
                                                          SE = round(sd(MSE_that_iteration, na.rm=T)/sqrt(simIterations), 3))

performanceDisplayName <- twoEntriesToValParenth(names(performanceTable)[2:3])

performanceTable[, performanceDisplayName] <- apply(performanceTable[,2:dim(performanceTable)[2]], MARGIN=1, 
                                              FUN=function(x) twoEntriesToValParenth(x))

performanceForDisplay <- performanceTable[, c(1,4)]

### Table of performance metrics
performanceTableLatex <- kable(performanceForDisplay, "latex", booktabs = T)
                                                    
fixed_table <- str_replace_all(performanceTableLatex, pattern = "textbackslash\\{\\}",
                               replace='')
      
fixed_table <- str_replace_all(fixed_table, pattern = '\\\\(?=(\\{|\\}))', replace='') 

fixed_table <- str_replace_all(fixed_table, pattern = "toprule|midrule|bottomrule",
                               replace='hline')                                                    
                                                    

### Boxplot of performance values
performancePlot <- (ggplot(performanceData, aes(x=Method, y=MSE_that_iteration, fill=Method)) + 
                      geom_boxplot() + theme_bw(base_size=20) + 
                      theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank()) +ylab("MSE")
                    + ggtitle("MSE for Score Function Estimation"))

rootOutputDir <- "../reports/score_performance_reports"

if(!dir.exists(rootOutputDir)){
  
  dir.create(rootOutputDir, recursive=TRUE)
  
  
}

                                                    
                                                    
### THIS SECTION WILL CREATE A FILE DIRECTORY FOR SAVING SIMULATION DATA,
### CREATED ONE DIRECTORY ABOVE LOCATION OF THIS FILE. 

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


ggsave(plot=performancePlot, filename=paste(outputDirectory, "Performance Boxplots.pdf", sep='/'), 
       units='in', width=8, height=8)


ggsave(plot=distBoxPlots, filename=paste(outputDirectory, "Unmixing Matrix Distance Boxplots.pdf", sep='/'), 
       units='in', width=8, height=8)
                                                    
ggsave(plot=coordinateErrorExample, filename=paste(outputDirectory, "FastICA vs us Error Example.pdf", sep='/'), 
       units='in', width=8, height=8)
                                                                                         
ggsave(plot=coordinateErrorExample , 
       filename=paste(outputDirectory, "Coordinate Error Histogram Example.pdf", sep='/'), 
       units='in', width=8, height=8)


ggsave(plot=coordinateDistributionExample , 
       filename=paste(outputDirectory, "Unmixed Coordinate Distribution Example.pdf", sep='/'), 
       units='in', width=8, height=8)
                                                    
ggsave(plot=unmixedCoordinateExamplePlot , 
       filename=paste(outputDirectory, "Unmixed Score Function Coordinate Example.pdf", sep='/'), 
       units='in', width=8, height=8)                


fileConn<-file(paste(outputDirectory, "Performance MSE Table.txt", sep='/'))
writeLines(fixed_table, fileConn)
close(fileConn)

fileConn<-file(paste(outputDirectory, "Unmixing Matrix Distance Table.txt", sep='/'))
writeLines(distanceTableLatex, fileConn)
close(fileConn)

write.csv(sim_params, file=paste(outputDirectory, 'parameters.csv', sep='/'))

