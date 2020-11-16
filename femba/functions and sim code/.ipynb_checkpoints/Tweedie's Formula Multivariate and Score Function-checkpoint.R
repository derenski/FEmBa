####### This defines the function the user interacts with directly to perform
####### bias correction.
library(caret)
library(splines2)
library(matrixcalc)
library(pracma)
library(expm)
library(doParallel)
library(foreach)

source("Function Data Simulation Parameter Makers.R")


#####################################################
#########################################################


##################### 


tweedieCorrectionNonJointICA <- function(X, Sigma_gamma,  ### Tweedie without joint unmixing and score function estimation
                                     functionalBasis,  
                                     transformation="none",
                                     gridStep = .001, 
                                     lambdaGrid=10^seq(-6, 6, 1),
                                     numberOfFolds = 5, 
                                     numberOfKnotsScoreFunction = 8){
    
    Sigma_gamma = .5*(Sigma_gamma+t(Sigma_gamma))
    
    scoreListToVector <- function(aScoreList){
    
    ### The output. Its arguments are a numberic vector, and an argument for desired order of derivative
    scoreVector <- function(theta, derivativeOrder=0){
      
      if (!is.matrix(theta)){
        
        theta <- array(theta, dim=c(length(theta), 1))
        
      }
      
      t(sapply(1:dim(theta)[1], FUN=function(x) aScoreList[[x]](theta[x,], derivativeOrder=derivativeOrder)))
      
    }
    
    return(scoreVector)
    
    
  }
  
  if (any(is.na(X))){
    
    stop('There are missing values in X.')

  }
  
  if (any(is.na(Sigma_gamma))){
    
    stop('There are missing values in the variance covariance matrix.')
    
  }
  
  p <- dim(X)[1]
  
  n <- dim(X)[2]
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis)/dim(functionalBasis)[1])
  
  sigmaTilde <- Sigma_gamma %*% basisCovarianceMatrix %*% Sigma_gamma
               
  sigmaTildeMinusHalf <- matrixToPower(sigmaTilde, -.5)
               
  covThetaMinusHalf <- matrixToPower(cov(t(X)), -.5)
               
  XSigmaTildeCov <- matrixToPower(sigmaTilde, .5) %*% covThetaMinusHalf %*% X

  if (transformation=="ica"){
      
    fastICAInfo <- fastICA(t(XSigmaTildeCov), n.comp = p)
      
    WThetaTildeEst <- t(fastICAInfo$K %*% fastICAInfo$W)
        
    UEst <- WThetaTildeEst %*% solve(sigmaTildeMinusHalf)
     
    WEst <-   UEst %*% covThetaMinusHalf
        
    Z <- WThetaTildeEst %*% XSigmaTildeCov
      
      
  }else if (transformation=='decorrelate'){
        
    WThetaTildeEst <- sigmaTildeMinusHalf 
     
    WEst <- covThetaMinusHalf
        
    Z <- WThetaTildeEst %*% XSigmaTildeCov
    
  }else if (transformation=="none"){
      
    WThetaTildeEst <- sigmaTildeMinusHalf
      
    WEst <- diag(1/diag(cov(t(X))))
    
    Z <- WThetaTildeEst %*% XSigmaTildeCov

  }else{
      
      stop("Invalid transformation, choose between one of (ica, decorrelate, none).")
      
  }  
               
  oneVec <- rep(1, dim(WThetaTildeEst)[1])
  
  cWs <- apply(WThetaTildeEst, MARGIN=1, 
               FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
               

  psi <- sapply(seq(1, dim(Z)[1]), FUN=function(i) scoreFunctionUpdate(Z=Z[i,], 
                              cW=cWs[i], lambdaGrid=lambdaGrid, numKnots=numberOfKnotsScoreFunction,
                              gridStep=gridStep, nFolds=numberOfFolds), simplify = F)
    
  ### Updated score vector
  psiVector <- scoreListToVector(aScoreList=psi)
  
  tweedieCorrection <- X + Sigma_gamma %*% t(WEst) %*% psiVector(Z)
  
  outputList <- list(tweedieEstimates=tweedieCorrection, W=WEst)
 
  return(outputList)
  
}






































tweediesFormulaOracle <- function(X, W, U, Sigma_gamma, functionalBasis, gridStep = .001, 
                                                  lambdaGrid=10^seq(-6, 6, 1),
                                                  numberOfFolds = 5, 
                                                  numberOfKnotsScoreFunction = 8,
                                                  maxIterations=100){
  
  scoreListToVector <- function(aScoreList){
    
    ### The output. Its arguments are a numberic vector, and an argument for desired order of derivative
    scoreVector <- function(theta, derivativeOrder=0){
      
      if (!is.matrix(theta)){
        
        theta <- array(theta, dim=c(length(theta), 1))
        
      }
      
      t(sapply(1:dim(theta)[1], FUN=function(x) aScoreList[[x]](theta[x,], derivativeOrder=derivativeOrder)))
      
    }
    
    return(scoreVector)
    
    
  }
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis)/dim(functionalBasis)[1])
  
  sigmaTilde <- Sigma_gamma %*% basisCovarianceMatrix %*% Sigma_gamma
  
  XSigmaTildeCov <- matrixToPower(sigmaTilde, .5) %*% matrixToPower(cov(t(X)), -.5) %*% X
  
  WthetaTilde <- U %*% matrixToPower(sigmaTilde, -.5)
  
  oneVec <- array(rep(1, dim(sigmaTilde)[1]), dim=c(dim(sigmaTilde)[1], 1))
  
  cWs <- apply(WthetaTilde, MARGIN=1, 
               FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
  
  XUnmixed <- WthetaTilde %*% XSigmaTildeCov
  
  trueScoreCoordinates <- sapply(seq(1, dim(XUnmixed)[1]), FUN=function(i) scoreFunctionUpdate(Z=XUnmixed[i,], 
                cW=cWs[i], lambdaGrid=lambdaGrid, numKnots=numberOfKnotsScoreFunction, 
                gridStep=gridStep, nFolds=numberOfFolds), simplify = F)
  
  TrueUnmixedScoreVector <- scoreListToVector(aScoreList=trueScoreCoordinates)
  
  tweedieCorrectedXs <- X + Sigma_gamma %*% t(W) %*% TrueUnmixedScoreVector(XUnmixed)
  
  outputInfo <- list(tweedieEstimates=tweedieCorrectedXs, W=W)                        
                                 
  return(outputInfo)
  
}


















matrixToPower <- function(aMatrix, power){ ### Utility for calculating matrices to given powers
                                           ### Assumes square, symmetric matrix
  matrixEigenDecomp <- eigen(aMatrix)
  
  matrixPowered <- matrixEigenDecomp$vectors %*% diag(matrixEigenDecomp$values^power) %*% t(matrixEigenDecomp$vectors)
  
  return(matrixPowered)
  
}



tweedieCorrectionWithICARandomRestart <- function(X, Sigma_gamma,  ### Tweedie's formula with estimtated score function.
                                     functionalBasis,        ### Score function estimated with ICA approach
                                     gridStep = .001, 
                                     lambdaGrid=10^seq(-6, 6, 1),
                                     numberOfFolds = 5, 
                                     numberOfKnotsScoreFunction = 8,
                                     numRestarts=3,
                                     numFoldsRestarts=2,
                                  algorithmSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06),
                                  updateUSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06)){
    
    if (numRestarts > 1 & numFoldsRestarts==1){
        
        stop("Must Have Validation Set (numFoldsRestarts > 1).")
        
    }
    
    options(warn=-1)
    
  evaluateRisk <- function(thetas, scoreFunction, SigmaTilde, W){
      
      uncoupoledScoreValues <- scoreFunction(W %*% thetas, derivativeOrder=0)
      
      uncoupledScoreGradients <- scoreFunction(W %*% thetas, derivativeOrder=1)
      
      WSigmaWT <- W %*% SigmaTilde %*% t(W)
      
      oneVec <- matrix(1, nrow=dim(W)[2], ncol=1)
      
      cWs <- apply(W, MARGIN=1, 
               FUN = function(x) t(oneVec) %*% (SigmaTilde * (x %*% t(x))) %*% oneVec)
      
      firstTerms <- mean(apply(uncoupoledScoreValues, MARGIN=2, FUN = function(x) t(x) %*% WSigmaWT %*% x))
                          
      secondTerms <- as.numeric(t(cWs) %*% matrix(rowMeans(uncoupledScoreGradients)))
                          
      riskValue <- firstTerms+2*secondTerms
                               
      return(riskValue)
      
  }  
    
  restartValidationFolds <- createFolds(1:dim(X)[2], k =numFoldsRestarts) 
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis
                             /dim(functionalBasis)[1]))
  
  sigmaTilde <- (Sigma_gamma %*% 
                   basisCovarianceMatrix %*% Sigma_gamma) ### Sigma tilde
  
  covTheta <- cov(t(X)) ### Estimated covariance of data
  
  ### Covariance of data to the -1/2 power
  covThetaMinusHalf <- matrixToPower(covTheta, -.5)
  
  ### Sigma tilde to the -1/2 power
  sigmaTildeMinusHalf <- matrixToPower(sigmaTilde, -.5)                            
  
  ### Initialization of U, where U is in so(p) 
                               
  UInits <- lapply(1:numRestarts, FUN=function(x) gramSchmidt(array(rnorm(dim(X)[1]^2), dim=rep(dim(X)[1], 2)))$Q)   

  risksRestarts <- foreach(UInit=UInits, .combine=c, .packages=c("doParallel", "foreach", "splines2", "caret"),
                          .export=c("scoreFunctionWithICA", "matrixToPower", "scoreFunctionUpdate",
                                   "unmixingMatrixUpdate")) %dopar% {                             
      
      risksThisRestart <- rep(NA, length(restartValidationFolds)  )    

      for (restartFold in restartValidationFolds){  

          X_train <- X[, restartFold]

          X_validation <- X[, -restartFold]

          ### Estimates the unmixing matrix, and the associated score functions
          scoreFunctionAndUnmixingMatrix <- scoreFunctionWithICA(thetaMatrix=X_train, sigmaTilde=sigmaTilde, 
                                          UInit=UInit, lambdaGrid=lambdaGrid, 
                                          numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                          nFolds=numberOfFolds, maxIterations=maxIterations,
                                          algorithmSpecs=algorithmSpecs, updateUSpecs=updateUSpecs)

          ### The U in so(p) that is used in calculating W, and W tilde
          finalU <- scoreFunctionAndUnmixingMatrix$U

          finalW <- finalU %*% covThetaMinusHalf

          ### The score function vector
          scoreVectors <- scoreFunctionAndUnmixingMatrix$score_vector
          
          if (numFoldsRestarts > 1){

              riskThisFold <- evaluateRisk(thetas=X_validation, scoreFunction=scoreVectors, SigmaTilde=sigmaTilde,                 W=finalW)
              
          }else{
                 
              riskThisFold <- -1*Inf
              
          }
          
          risksThisRestart <- c(risksThisRestart, riskThisFold)

          }
      
          mean(risksThisRestart, na.rm=T)
      
    }
                               
    bestInitialization <- UInits[[which.min(risksRestarts)]]
                               
    scoreFunctionAndUnmixingMatrix <- scoreFunctionWithICA(thetaMatrix=X, sigmaTilde=sigmaTilde, 
                                          UInit=bestInitialization, lambdaGrid=lambdaGrid, 
                                          numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                          nFolds=numberOfFolds, maxIterations=maxIterations,
                                          algorithmSpecs=algorithmSpecs, updateUSpecs=updateUSpecs)

          ### The U in so(p) that is used in calculating W, and W tilde
    finalU <- scoreFunctionAndUnmixingMatrix$U

    finalW <- finalU %*% covThetaMinusHalf

          ### The score function vector
    scoreVectors <- scoreFunctionAndUnmixingMatrix$score_vector
                               
          ### Uncoupled data based on our estimated unmixing matrix
    uncoupledThetaTildes <- finalW %*% X

          ### Values for score function of unmixed coordinates
    scoreValues <- t(finalW) %*% scoreVectors(uncoupledThetaTildes)

          ### The tweedie correction
    tweedieCorrection <- X + Sigma_gamma %*% scoreValues       
                   
    outputInfo <- list(tweedieEstimates=tweedieCorrection, W=finalW)
                               
    return(outputInfo)
  
}
                   
                   
###################### With warm start provided by FastICA
                   
tweedieCorrectionWithICAWarmStart <- function(X, Sigma_gamma,  ### Tweedie's formula with estimtated score function.
                                     functionalBasis,        ### Score function estimated with ICA approach
                                     gridStep = .001, 
                                     lambdaGrid=10^seq(-6, 6, 1),
                                     numberOfFolds = 5, 
                                     numberOfKnotsScoreFunction = 8,
                                  algorithmSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06),
                                  updateUSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06)){
    
    options(warn=-1)   
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis
                             /dim(functionalBasis)[1]))
  
  sigmaTilde <- (Sigma_gamma %*% 
                   basisCovarianceMatrix %*% Sigma_gamma) ### Sigma tilde
  
  covTheta <- cov(t(X)) ### Estimated covariance of data
  
  ### Covariance of data to the -1/2 power
  covThetaMinusHalf <- matrixToPower(covTheta, -.5)
  
  ### Sigma tilde to the -1/2 power
  sigmaTildeMinusHalf <- matrixToPower(sigmaTilde, -.5)                            
  
  ### Initialization of U, where U is in so(p) 
   XSigmaTildeCov <- matrixToPower(sigmaTilde, .5) %*% covThetaMinusHalf %*% X

   fastICAInfo <- fastICA(t(XSigmaTildeCov), n.comp = dim(Sigma_gamma)[1])
    
   WThetaTildeEst <- t(fastICAInfo$K %*% fastICAInfo$W)
    
   UWarmStart <- WThetaTildeEst %*% solve(sigmaTildeMinusHalf)
    
   # XUncorr <- covThetaMinusHalf %*% X
    
   # fastICAInfo <- fastICA(t(XUncorr), n.comp = dim(Sigma_gamma)[1])
    
  # UWarmStart <- t(fastICAInfo$W)

  scoreFunctionAndUnmixingMatrix <- scoreFunctionWithICA(thetaMatrix=X, sigmaTilde=sigmaTilde, 
                                          UInit=UWarmStart, lambdaGrid=lambdaGrid, 
                                          numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                          nFolds=numberOfFolds, maxIterations=maxIterations,
                                          algorithmSpecs=algorithmSpecs, updateUSpecs=updateUSpecs)

          ### The U in so(p) that is used in calculating W, and W tilde
    finalU <- scoreFunctionAndUnmixingMatrix$U

    finalW <- finalU %*% covThetaMinusHalf
    
    finalWTilde <- finalU %*% sigmaTildeMinusHalf

          ### The score function vector
    scoreVectors <- scoreFunctionAndUnmixingMatrix$score_vector
                               
          ### Uncoupled data based on our estimated unmixing matrix
    #uncoupledThetaTildes <- finalW %*% X
    
    uncoupledThetaTildes <- finalWTilde %*% XSigmaTildeCov

          ### The tweedie correction
    tweedieCorrection <- X + Sigma_gamma %*% t(finalW) %*% scoreVectors(uncoupledThetaTildes)
    
    
    outputInfo <- list(tweedieEstimates=tweedieCorrection, W=finalW)
                               
    return(outputInfo)
  
}
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   


### Names of algorithmSpecs, updateUSpecs: 
## c("min_iterations", "max_iterations", "min_distance")

### Function that runs score function estimation algorithm
scoreFunctionWithICA <- function(thetaMatrix, sigmaTilde, UInit, lambdaGrid=10^seq(-6, 6, 1), numKnots=5, 
            gridStep=.01, nFolds=5, maxIterations=200, maxIterationsUpdateU=100,
            algorithmSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06),
            updateUSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06)){ ## thetaMatrix is p times n
    
        ## You can use this to evaluate the risk for the UNTRANSFORMED data (compare fastICA's risk to yours)
    evaluateRisk <- function(thetas, unmixedScoreFunction, SigmaTilde, W){
      
      uncoupoledScoreValues <- unmixedScoreFunction(W %*% thetas, derivativeOrder=0)
      
      uncoupledScoreGradients <- unmixedScoreFunction(W %*% thetas, derivativeOrder=1)
      
      WSigmaWT <- W %*% SigmaTilde %*% t(W)
      
      oneVec <- matrix(1, nrow=dim(W)[2], ncol=1)
      
      cWs <- apply(W, MARGIN=1, 
               FUN = function(x) t(oneVec) %*% (SigmaTilde * (x %*% t(x))) %*% oneVec)
      
      firstTerms <- mean(apply(uncoupoledScoreValues, MARGIN=2, FUN = function(x) t(x) %*% WSigmaWT %*% x))
                          
      secondTerms <- as.numeric(t(cWs) %*% matrix(rowMeans(uncoupledScoreGradients)))
                          
      riskValue <- firstTerms+2*secondTerms
                               
      return(riskValue)
      
  }
   
  ### Converts list of functions to a vector function.
  scoreListToVector <- function(aScoreList){
    
    ### The output. Its arguments are a numberic vector, and an argument for desired order of derivative
    scoreVector <- function(theta, derivativeOrder=0){
      
      if (!is.matrix(theta)){
        
        theta <- array(theta, dim=c(length(theta), 1))
        
      }
      
      t(sapply(1:dim(theta)[1], FUN=function(x) aScoreList[[x]](theta[x,],
                                                                derivativeOrder=derivativeOrder)))
      
    }
    
    return(scoreVector)
    
    
  }
  
  ### Covariance of Data
  covTheta <- cov(t(thetaMatrix))
  
  ### Covariance of Data to the -1/2 power
  covThetaMinusHalf <- matrixToPower(covTheta, -.5)
  
  ### Sigma tilde to the -1/2 power
  sigmaTildeMinusHalf <- matrixToPower(sigmaTilde, -.5)
  
  ### Rotation of data so it has covariance Sigma tilde
  thetasSigmaTildeCov <- solve(sigmaTildeMinusHalf) %*% covThetaMinusHalf %*% thetaMatrix
  
  oneVec <- matrix(1, nrow=dim(UInit)[2], ncol=1)
  
  ### Initial W tilde
  WtThetaTilde <- UInit %*% sigmaTildeMinusHalf
  
  ### Intial W
  Wt <- UInit %*% covThetaMinusHalf
  
  ### Constants needed for estimation score vector
  cWs <- apply(WtThetaTilde, MARGIN=1, 
               FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
  
  ### Initialization of uncoupled data
  Zt <- WtThetaTilde %*% thetasSigmaTildeCov
  
  ### Initialization of score functions
  psi_t <- sapply(seq(1, dim(Zt)[1]), FUN=function(i) scoreFunctionUpdate(Z=Zt[i,], cW=cWs[i], 
        lambdaGrid=lambdaGrid, numKnots=numKnots, gridStep=gridStep, nFolds=nFolds), simplify = F)
  
  ### Convert list to vector function
  psi_tVector <- scoreListToVector(aScoreList=psi_t)
  
  ### Score function values for unmixed data
  unmixedScoreVals_t <- psi_tVector(Zt, derivativeOrder = 0)
  
  ### Score function derivative values for unmixed data
  unmixedScoreGradient_t <- psi_tVector(Zt, derivativeOrder = 1)
  
  ### Risk value at initialization
  risk_t <- evaluateRisk(thetas=thetaMatrix, unmixedScoreFunction=psi_tVector, SigmaTilde=sigmaTilde,
                        W=Wt)
  
  mixedScoreVals_t <- t(Wt) %*% unmixedScoreVals_t
  
  print(paste("Initial Risk:", risk_t))                
                  
  for (iterator in 1:algorithmSpecs["max_iterations"]){

    ### The prewhitening done for updating the densities and updating W are slightly different.
    ### This is because to update W we further decomepose it into two matrices
    
    ### Current U in so(p). We first update Ut and then update the score functions based on Ut.
    Ut <- WtThetaTilde %*% solve(sigmaTildeMinusHalf)
    
    ### The updated Ut
    UtPlusOne <- unmixingMatrixUpdate(thetaStarMatrix= thetaMatrix, sigmaTilde=sigmaTilde,
                                      UInit=Ut, 
                                    scoreFunctionVectorInit=psi_tVector, 
                                    minIterationsU = updateUSpecs["min_iterations"],
                                    maxIterationsU = updateUSpecs["max_iterations"],
                                    minChangeU = updateUSpecs["min_change"])
    
    ### The new W wilde
    WtPlusOneThetaTilde <- UtPlusOne %*% sigmaTildeMinusHalf
    
    ### The new W
    WtPlusOne <- UtPlusOne %*% covThetaMinusHalf
    
    ### Constants needed for score function estimation
    cWs <- apply(WtPlusOneThetaTilde, MARGIN=1, 
                 FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
    
    ### Unmixed data based on new W tilde
    ZtPlusOne <- WtPlusOneThetaTilde %*% thetasSigmaTildeCov
    
    ### Updated score functions
    psi_tPlusOne <- sapply(seq(1, dim(ZtPlusOne)[1]), FUN=function(i) scoreFunctionUpdate(Z=ZtPlusOne[i,], 
                              cW=cWs[i], lambdaGrid=lambdaGrid, numKnots=numKnots,
                              gridStep=gridStep, nFolds=nFolds), simplify = F)
    
    ### Updated score vector
    psi_tPlusOneVector <- scoreListToVector(aScoreList=psi_tPlusOne)
    
    ### Values needed to calculate risk
    unmixedScoreVals_tPlusOne <- psi_tPlusOneVector(ZtPlusOne, derivativeOrder = 0)
    
    unmixedScoreGradient_tPlusOne <- psi_tPlusOneVector(ZtPlusOne, derivativeOrder = 1)
    
    ### The risk at the updated unmixing matrix and score vector. Calculated based on the unmixed score vector.
    risk_tPlusOne <- evaluateRisk(thetas=thetaMatrix, unmixedScoreFunction=psi_tPlusOneVector,
                                  SigmaTilde=sigmaTilde, W=WtPlusOne)
    
    mixedScoreVals_tPlusOne <- t(WtPlusOne) %*% unmixedScoreVals_tPlusOne

    ### The absolute percentage change in the risk
    Change <- abs(risk_t-risk_tPlusOne)/abs(risk_t)
    
    overMinIterationsU <- iterator > algorithmSpecs["min_iterations"]
    
    closerThanDistU <- Change < algorithmSpecs["min_change"]
    
    if (overMinIterationsU & closerThanDistU){
      
      break
    
    }else{
      
      ### Updating the needed parameters for the next iteration
      
      WtThetaTilde <- WtPlusOneThetaTilde
      
      Wt <- WtPlusOne
      
      Zt <- ZtPlusOne
      
      psi_t <- psi_tPlusOne
      
      psi_tVector <- psi_tPlusOneVector
      
      unmixedScoreVals_t <- unmixedScoreVals_tPlusOne
      
      mixedScoreVals_t <- mixedScoreVals_tPlusOne
      
      risk_t <- risk_tPlusOne
      
    }
    
    #print(Change)
    
  }
    
  ### Return final U needed for calculating unmixing matrices, and the score vector
                           
    print(paste("Final Risk:", risk_tPlusOne))                       
    
    UandScoreVec <- list(UtPlusOne, psi_tPlusOneVector)
    
    names(UandScoreVec) <- c("U", "score_vector")

  return(UandScoreVec)
  
}
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           


### Update Unmixing Matrix

unmixingMatrixUpdate <- function(thetaStarMatrix, UInit, sigmaTilde,
                                 scoreFunctionVectorInit, minIterationsU=10, maxIterationsU=100,
                                 minChangeU=1e-06){
    
  Gamma <- cov(t(thetaStarMatrix))
    
  thetaTildeMatrix <- matrixToPower(Gamma, -.5) %*% thetaStarMatrix
  
  ### Note that A is a symmetric matrix
  A <- matrixToPower(Gamma, -.5) %*% sigmaTilde %*% matrixToPower(Gamma, -.5)
    
  ### The value of the risk
  rHat <- function(uRisk, ourScoreVector, thetaTildeMatrix, A){ 
      
    UtAUrisk <- uRisk %*% A %*% t(uRisk)
    
    uncoupledThetaTildesUrisk <- uRisk %*% thetaTildeMatrix
      
    riskScoreVals <- ourScoreVector( uncoupledThetaTildesUrisk, derivativeOrder=0)
      
    riskScoreValsPrime <- ourScoreVector( uncoupledThetaTildesUrisk, derivativeOrder=1)
    
    firstPart <- sum(t(riskScoreVals) %*% UtAUrisk %*% riskScoreVals)/dim(riskScoreVals)[2]
      
    secondPart <- mean(2*array(diag(UtAUrisk), dim=c(1, dim(riskScoreValsPrime)[1])) %*% riskScoreValsPrime)
    
    return(firstPart+secondPart)
    
  }
  
  ### The gradient of the risk
  nablaRHat <- function(uRisk, thetaTildeMatrix, scoreFunctionVector, A){
      
    crossTermCalculatorForU <- function(scoreVector, theU){
        
        reverseIdentityMat <- array(1, dim=rep(length(scoreVector), 2))-diag(rep(1, 
                            length(scoreVector)))
        
        scoreVectorMadeMatrixWithMissingDiag <- sapply(scoreVector, FUN=rep,
                                                       length(scoreVector))*reverseIdentityMat
        
        neededTerm <- scoreVectorMadeMatrixWithMissingDiag %*% theU
        
        scoreVectorAsMatrix <- t(sapply(scoreVector, FUN=rep, dim(uRisk)[1]))
        
        return(scoreVectorAsMatrix*neededTerm)
        
    }    
      
    crossTermCalculatorForThetaTilde <- function(scoreVector, scoreVectorPrime, UtAU){
        
        reverseIdentityMat <- array(1, dim=rep(length(scoreVector), 2))-diag(rep(1,
                                    length(scoreVector)))

        scoreVectorMadeMatrixWithMissingDiag <- UtAU*sapply(scoreVector, FUN=rep,
                                                       length(scoreVector))*reverseIdentityMat
        
        scoreVectorPrimeAsMatrix <- t(sapply(scoreVectorPrime, FUN=rep, dim(uRisk)[1]))
        
        return(2*scoreVectorPrimeAsMatrix*scoreVectorMadeMatrixWithMissingDiag)
        
    }     
      
      
    
    uncoupledThetaTildesUrisk <- uRisk %*% thetaTildeMatrix
      
    scoreVals <- scoreFunctionVector(uncoupledThetaTildesUrisk, derivativeOrder=0)
      
    scoreValsPrime <- scoreFunctionVector(uncoupledThetaTildesUrisk, derivativeOrder=1)
      
    scoreValsDoublePrime <- scoreFunctionVector(uncoupledThetaTildesUrisk, derivativeOrder=2)
     
    UtAU <- uRisk %*% A %*% t(uRisk)
      
      
    ### Terms that get multiplied to U, no cross terms
    UkTermsNonMatrix <- rowMeans(2*(scoreVals^2)+4*scoreValsPrime)
      
    UkTermsMatrix <- t(sapply(UkTermsNonMatrix, FUN=rep, dim(uRisk)[1]))
      
    firstPart <- UkTermsMatrix * (uRisk %*% t(A)) ### Need row means
      
    ### Terms that get multiplied to theta tilde  
    thetaTildeTermsMatrix <- t(sapply(diag(UtAU), 
                                      FUN=rep, 
                                      dim(scoreVals)[2]))*(2*scoreVals*scoreValsPrime+
                                                                    2*scoreValsDoublePrime)  
      
    secondPart <- (thetaTildeTermsMatrix %*% t(thetaTildeMatrix))/dim(thetaTildeMatrix)[2]
      
      
    ### Terms that get multiplied to U, with cross terms  
    thirdPartUnAggregatedU <- apply(scoreVals, MARGIN=2, FUN= crossTermCalculatorForU, 
                                   theU=uRisk %*% t(A))
      
    thirdPartU <- Reduce("+", thirdPartUnAggregatedU) / length(thirdPartUnAggregatedU)
      
      
    thirdPartUnAggregatedTheta <- sapply(1:dim(scoreVals)[2], FUN=function(x)
        crossTermCalculatorForThetaTilde(scoreVector=scoreVals[,x],
                                        scoreVectorPrime=scoreValsPrime[,x],
                                        UtAU=UtAU))
                                         
    thetasExpanded <- sapply(1:dim(thetaTildeMatrix)[2], 
                             FUN=function(x) sapply(x, FUN=rep, length(x)))
                             
    listOfThetaCrossTerms <- sapply(1:length(thetasExpanded), FUN=function(x)
      thirdPartUnAggregatedTheta[[x]] %*% thetasExpanded[[x]]  )
                                         
    thirdPartThetaCrossTerms <- Reduce("+", listOfThetaCrossTerms) / length(listOfThetaCrossTerms)
      
      
    nablaValues <- firstPart + secondPart + thirdPartU + thirdPartThetaCrossTerms

    return(nablaValues)
    
  } 
  
  scoreFunctionVector <- scoreFunctionVectorInit
  
  UtMinus1 <- UInit
  
  dist <- 1
  
  alpha <- .2
  
  etaGrid <-10*c(.5^seq(0, 15, 1), 0)
  
  thingy <- 1
  
  for (iteratorU in 1:maxIterationsU){
    
    ### Unmix decorrelated data
    uncoupledThetaTildes <- UtMinus1 %*% thetaTildeMatrix 
    
    ### Gradient of risk
    nablaRHat_t <- nablaRHat(uRisk=UtMinus1, thetaTildeMatrix=thetaTildeMatrix, scoreFunctionVector,A=A)
    
    ### Geodesic gradient
    nablaRHatTilde <- nablaRHat_t %*% t(UtMinus1) - UtMinus1 %*% t(nablaRHat_t)
    
    Ht <- nablaRHatTilde
    
    ### We normalize the gradient as described in the notes
    Ht <- sqrt(2)*Ht/norm(Ht, "F")
    
    expEtaHs <- sapply(etaGrid, FUN=function(x) expm::expm(-1*x * Ht, order=16), simplify = F)
    
    potentialUs <- lapply(expEtaHs, FUN=function(x) x %*% UtMinus1)
    
    riskSatisfactionTerms <- function(x, y){

      return(c(rHat(uRisk=x, ourScoreVector=scoreFunctionVector, thetaTildeMatrix=thetaTildeMatrix, A=A), 
               rHat(uRisk=UtMinus1, ourScoreVector=scoreFunctionVector,thetaTildeMatrix=thetaTildeMatrix, A=A)- 
               alpha*.5*y* (norm(Ht, "F")^2)))
      
    }
    
    neededRiskTerms <- mapply(riskSatisfactionTerms, x=potentialUs, y=etaGrid)
    
    doesItSatisfy <- apply(neededRiskTerms, 2, FUN = function(x) x[1] >= x[2])
    
    if(all(doesItSatisfy)){
      ### If we can't improve, we choose our step to be the identity matrix (so no change)
      expEtaH <- expEtaHs[[length(expEtaHs)]]
      
    }else{
      ### Choose the risk minimizer 
      expEtaH <- expEtaHs[[which.min(neededRiskTerms[1,])]]
    
    }
    ### Update U
    Ut <- expEtaH %*% UtMinus1
    
    ### Percentage change in risk
    dist <- (rHat(uRisk=UtMinus1, ourScoreVector=scoreFunctionVector,thetaTildeMatrix=thetaTildeMatrix,A=A)-
                   rHat(uRisk=Ut, ourScoreVector=scoreFunctionVector,  
                thetaTildeMatrix, A=A))/abs(rHat(uRisk=UtMinus1, ourScoreVector=scoreFunctionVector,
                                            thetaTildeMatrix=thetaTildeMatrix, A=A))

    overMinIterationsU <- iteratorU > minIterationsU
    
    overMaxIterationsU <- iteratorU > maxIterationsU
    
    closerThanDistU <- dist < minChangeU
      
    if (overMinIterationsU & closerThanDistU){
        
      break
        
    }else{ 
        
      UtMinus1 <- Ut
        
      }
      
  }
  
  return(Ut)
  
}





                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           


#### Update Score function

scoreFunctionUpdate <- function(Z, cW, lambdaGrid=10^seq(-6, 6, 1), numKnots=10, gridStep=.01,
                                    nFolds=5){
  
  ### Calculates the risk
  riskCalculator <- function(Z, betaHat, basis, lambda, cW, gridStep){
    
    n <- length(Z)
    
    sOfZ <- predict(basis, newx = Z)
    
    sMatrix <- base::crossprod(sOfZ)/n
    
    basisFirstDerivative <- deriv(basis, derivs = 1L)
    
    sPrimeOfZ <- predict(basisFirstDerivative, newx = Z)
    
    sPrimeVec <- matrix(rowMeans(t(sPrimeOfZ)), nrow=dim(sPrimeOfZ)[2], ncol=1)
    
    basisSecondDerivative <- deriv(basis, derivs = 2L)
    
    omega <- gridStep*base::crossprod(basisSecondDerivative)
    
    riskValue <- t(betaHat) %*% (sMatrix+lambda*omega) %*% betaHat + 2*cW * t(sPrimeVec) %*% betaHat
    
    return(riskValue)
    
  }
  
  ### Calculates the basis coefficient vector for the score function
  betaCalculator <- function(Z, basis, lambda, cW, gridStep){
    
    n <- length(Z)
    
    sOfZ <- predict(basis, newx = Z)
    
    sMatrix <- base::crossprod(sOfZ)/n
    
    basisFirstDerivative <- deriv(basis, derivs = 1L)
    
    sPrimeOfZ <- predict(basisFirstDerivative, newx = Z)
    
    sPrimeVec <- matrix(rowMeans(t(sPrimeOfZ)), nrow=dim(sPrimeOfZ)[2], ncol=1)
    
    basisSecondDerivative <- deriv(basis, derivs = 2L)
    
    omega <- gridStep*base::crossprod(basisSecondDerivative)
    
    betaHat <- -cW*solve(sMatrix+lambda*omega) %*% sPrimeVec
    
    return(betaHat)
    
  }
  
  ### If bigger than 1, perform cross validation, otherwise estimate from given turning parameter
  CV <- length(lambdaGrid) > 1
    
    
  if (CV & (nFolds==1)){
      
      stop("Must have multiple folds if doing CV.")
      
  }
  
  ### Number of curves
  n <- length(Z) 
  
  num_knots_plus_1 <- numKnots+1 ### Give the number of knots, plus 1
  
  theProbs <- seq(1:(num_knots_plus_1-1))/(num_knots_plus_1)
  
  theKnots <- quantile(Z, probs = theProbs) ### Knots based on quantiles
  
  ### Making sure knots are not placed are extreme points in the data
  if (min(theKnots) == min(Z)){
    
    theKnots[which(theKnots==min(Z))] <- min(Z[which(Z != min(Z))])
    
    
  }
  
  if (max(theKnots) == max(Z)){
    
    theKnots[which(theKnots==max(Z))] <- min(Z[which(Z != max(Z))])
    
    
  }
  
  ### Generating the basis to be used to represent the score function
  theKnots <- unique(theKnots)
  
  zMin <- floor(min(Z, na.rm = T))
  
  zMax <- ceiling(max(Z, na.rm = T))
  
  zGrid <- seq(zMin, zMax,  gridStep)

  basis<- bSpline(seq(zMin,zMax,gridStep), degree = min(3,n), 
                         knots = theKnots, intercept = T) ### use knots
    
    #### number of knots = df-degree
   #### Specify knots to be at quantiles 
  
  df <- num_knots_plus_1+min(3,n)
  
  if (CV){ ### Cross-validation for tuning parameter estimation
    
    cvFolds <- createFolds(1:n, k =nFolds)
    
    ### Cross validation done in parallel
    cvErrorsForLambdas <- foreach(lamb=lambdaGrid,
                .packages=c('splines2', "caret"), .combine = 'c') %dopar% {
      
      risksThisLambda <- c()
      
      cvFolds <- createFolds(1:n, k =nFolds)
      
      for (fold in cvFolds){
      
        trainZ <- Z[-fold]
        
        valZ <- Z[fold]
      
        betaHatThisFoldThisLambda <- betaCalculator(Z=trainZ, basis=basis, lambda=lamb, cW=cW, gridStep=gridStep)
        
        riskThisFoldThisLambda <- riskCalculator(Z=valZ, betaHat=betaHatThisFoldThisLambda,
                                                 basis=basis, lambda=lamb, cW=cW, gridStep=gridStep)
        
        risksThisLambda <- c(risksThisLambda, riskThisFoldThisLambda)
      
      }
     
      mean(risksThisLambda)
       
    }
    
    chosenLambda <- lambdaGrid[which.min(cvErrorsForLambdas)]
    
  }else{
    
    chosenLambda <- lambdaGrid
    
  }
  
  betaHat <- betaCalculator(Z=Z, basis=basis, lambda=chosenLambda, cW=cW, gridStep=gridStep)
  
  ### The estimated score function that gets returned
  scoreFunction <- function(z, derivativeOrder=0){
    
    if (derivativeOrder==0){
      
      basisNow = basis
      
    }else{
    
    basisNow <- deriv(basis, derivs = derivativeOrder)
    
    }
    
    basisVec = predict(basisNow, newx = z)
    
    return(basisVec %*% betaHat)
    
  }
  
  force(scoreFunction)
  
  return(scoreFunction)
  
}












