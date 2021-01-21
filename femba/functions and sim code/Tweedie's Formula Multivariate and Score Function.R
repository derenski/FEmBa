####### This defines the function the user interacts with directly to perform
####### bias correction.
library(caret)
library(splines2)
library(MASS)
library(lava)
library(matrixcalc)
library(pracma)
library(expm)
library(doParallel)
library(foreach)
library(plyr)
library(fastICA)

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

  covThetaMinusHalf <- solve(sqrtm(cov(t(X))))

  if (transformation=="ica"){
      
    fastICAInfo <- fastICA(t(X), n.comp = p)
      
    WEst <- t(fastICAInfo$K %*% fastICAInfo$W)
        
    UEst <- WEst %*% solve(covThetaMinusHalf)
      
    Z <- WEst %*% X
      
  }else if (transformation=='decorrelate'){

    WEst <- covThetaMinusHalf
        
    Z <- WEst %*% X
    
  }else if (transformation=="none"){

    WEst <- diag(1/diag(cov(t(X))))
    
    Z <- WEst %*% X

  }else{
      
      stop("Invalid transformation, choose between one of (ica, decorrelate, none).")
      
  }  
               
  #oneVec <- rep(1, dim(WEst)[1])
  
  #cWs <- apply(WEst, MARGIN=1, 
  #             FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
               

  #psi <- sapply(seq(1, dim(Z)[1]), FUN=function(i) scoreFunctionUpdate(Z=Z[i,], 
  #                            cW=cWs[i], lambdaGrid=lambdaGrid, numKnots=numberOfKnotsScoreFunction,
  #                            gridStep=gridStep, nFolds=numberOfFolds), simplify = F)
    
  ### Updated score vector
  #psiVector <- scoreListToVector(aScoreList=psi)
               
               
  psiVector <- scoreFunctionUpdate(Z=Z, W=WEst, SigmaTilde=sigmaTilde, lambdaGrid=lambdaGrid, numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                    nFolds=numberOfFolds)
  
  tweedieCorrection <- X + Sigma_gamma %*% t(WEst) %*% psiVector(Z)
  
  outputList <- list(tweedieEstimates=tweedieCorrection, W=WEst)
 
  return(outputList)
  
}






































tweediesFormulaOracle <- function(X, W, U, Sigma_gamma, functionalBasis, gridStep = .001, 
                                                  lambdaGrid=10^seq(-6, 6, 1),
                                                  numberOfFolds = 5, 
                                                  numberOfKnotsScoreFunction = 8,
                                                  maxIterations=100){
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis)/dim(functionalBasis)[1])
  
  sigmaTilde <- Sigma_gamma %*% basisCovarianceMatrix %*% Sigma_gamma
  
  Wtheta <- W
  
  oneVec <- array(rep(1, dim(sigmaTilde)[1]), dim=c(dim(sigmaTilde)[1], 1))
  
  #cWs <- apply(Wtheta, MARGIN=1, 
  #             FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
  
  XUnmixed <- Wtheta %*% X
               
  TrueUnmixedScoreVector <- scoreFunctionUpdate(Z=XUnmixed, W=Wtheta, SigmaTilde=sigmaTilde, lambdaGrid=lambdaGrid, numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                    nFolds=numberOfFolds)
  
  tweedieCorrectedXs <- X + Sigma_gamma %*% t(W) %*% TrueUnmixedScoreVector(XUnmixed)
  
  outputInfo <- list(tweedieEstimates=tweedieCorrectedXs, W=W)                        
                                 
  return(outputInfo)
  
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
  covThetaMinusHalf <- solve(sqrtm(covTheta))                     
  
  ### Initialization of U, where U is in so(p) 
                               
  UInits <- lapply(1:numRestarts, FUN=function(x) gramSchmidt(array(rnorm(dim(X)[1]^2), dim=rep(dim(X)[1], 2)))$Q)   

  risksRestarts <- foreach(UInit=UInits, .combine=c, .packages=c("doParallel", "foreach", "splines2", "caret", "lava"),
                          .export=c("scoreFunctionWithICA", "scoreFunctionUpdate",
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
                                  updateUSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06),
                                             learningRateU=.2){
    
    options(warn=-1)   
  
  basisCovarianceMatrix <- ((t(functionalBasis) %*% ### Term calculated from basis used to model curves
                               functionalBasis
                             ))/dim(functionalBasis)[1]
  
  sigmaTilde <- (Sigma_gamma %*% 
                   basisCovarianceMatrix %*% Sigma_gamma) ### Sigma tilde
  
  covTheta <- cov(t(X)) ### Estimated covariance of data
  
  ### Covariance of data to the -1/2 power
  covThetaMinusHalf <- solve(sqrtm(covTheta))

   fastICAInfo <- fastICA(t(X), n.comp = dim(Sigma_gamma)[1])
    
   WEst <- t(fastICAInfo$K %*% fastICAInfo$W)
    
   UWarmStart <- WEst%*% solve(covThetaMinusHalf)
    
   # XUncorr <- covThetaMinusHalf %*% X
    
   # fastICAInfo <- fastICA(t(XUncorr), n.comp = dim(Sigma_gamma)[1])
    
  # UWarmStart <- t(fastICAInfo$W)

  scoreFunctionAndUnmixingMatrix <- scoreFunctionWithICA(thetaMatrix=X, sigmaTilde=sigmaTilde, 
                                          UInit=UWarmStart, lambdaGrid=lambdaGrid, 
                                          numKnots=numberOfKnotsScoreFunction, gridStep=gridStep,
                                          nFolds=numberOfFolds, maxIterations=maxIterations,
                                          algorithmSpecs=algorithmSpecs, updateUSpecs=updateUSpecs,
                                                        learningRateU=learningRateU)

          ### The U in so(p) that is used in calculating W, and W tilde
    finalU <- scoreFunctionAndUnmixingMatrix$U

    finalW <- finalU %*% covThetaMinusHalf

          ### The score function vector
    scoreVectors <- scoreFunctionAndUnmixingMatrix$score_vector
                               
          ### Uncoupled data based on our estimated unmixing matrix
    #uncoupledThetaTildes <- finalW %*% X
    
    uncoupledThetaTildes <- finalW %*% X

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
            updateUSpecs=c(min_iterations=10, max_iterations=100, min_change=1e-06),
                                learningRateU=.2){ ## thetaMatrix is p times n
    
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
  
  ### Covariance of Data
  covTheta <- cov(t(thetaMatrix))
  
  ### Covariance of Data to the -1/2 power
  covThetaMinusHalf <- solve(sqrtm(covTheta))
  
  oneVec <- matrix(1, nrow=dim(UInit)[2], ncol=1)
  
  ### Intial W
  Wt <- UInit %*% covThetaMinusHalf
  
  ### Constants needed for estimation score vector
  #cWs <- apply(Wt, MARGIN=1, 
  #             FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
  
  ### Initialization of uncoupled data
  Zt <- Wt %*% thetaMatrix
  
  ### Initialization of score functions
  #psi_t <- sapply(seq(1, dim(Zt)[1]), FUN=function(i) scoreFunctionUpdate(Z=Zt[i,], cW=cWs[i], 
  #      lambdaGrid=lambdaGrid, numKnots=numKnots, gridStep=gridStep, nFolds=nFolds), simplify = F)
  
  ### Convert list to vector function
  #psi_tVector <- scoreListToVector(aScoreList=psi_t)
                  
  psi_tVector <- scoreFunctionUpdate(Z=Zt, W=Wt, SigmaTilde=sigmaTilde, lambdaGrid=lambdaGrid, numKnots=numKnots, gridStep=gridStep,
                                    nFolds=nFolds)
  
  ### Score function values for unmixed data
  unmixedScoreVals_t <- psi_tVector(Zt, derivativeOrder = 0)
  
  ### Score function derivative values for unmixed data
  unmixedScoreGradient_t <- psi_tVector(Zt, derivativeOrder = 1)
  
  ### Risk value at initialization
  risk_t <- evaluateRisk(thetas=thetaMatrix, unmixedScoreFunction=psi_tVector, SigmaTilde=sigmaTilde,
                        W=Wt)
  
  mixedScoreVals_t <- t(Wt) %*% unmixedScoreVals_t              
                  
  for (iterator in 1:algorithmSpecs["max_iterations"]){

    ### The prewhitening done for updating the densities and updating W are slightly different.
    ### This is because to update W we further decomepose it into two matrices
    
    ### Current U in so(p). We first update Ut and then update the score functions based on Ut.
    Ut <- Wt %*% solve(covThetaMinusHalf)
    
    ### The updated Ut
    UtPlusOne <- unmixingMatrixUpdate(thetaStarMatrix= thetaMatrix, sigmaTilde=sigmaTilde,
                                      UInit=Ut, 
                                    scoreFunctionVectorInit=psi_tVector, 
                                    minIterationsU = updateUSpecs["min_iterations"],
                                    maxIterationsU = updateUSpecs["max_iterations"],
                                    minChangeU = updateUSpecs["min_change"],
                                     alpha=learningRateU)
    
    ### The new W
    WtPlusOne <- UtPlusOne %*% covThetaMinusHalf
    
    ### Constants needed for score function estimation
    #cWs <- apply(WtPlusOne, MARGIN=1, 
    #             FUN = function(x) t(oneVec) %*% (sigmaTilde * (x %*% t(x))) %*% oneVec)
    
    ### Unmixed data based on new W tilde
    ZtPlusOne <- WtPlusOne %*% thetaMatrix
    
    ### Updated score functions
   # psi_tPlusOne <- sapply(seq(1, dim(ZtPlusOne)[1]), FUN=function(i) scoreFunctionUpdate(Z=ZtPlusOne[i,], 
   #                           cW=cWs[i], lambdaGrid=lambdaGrid, numKnots=numKnots,
   #                           gridStep=gridStep, nFolds=nFolds), simplify = F)
    
    ### Updated score vector
    # psi_tPlusOneVector <- scoreListToVector(aScoreList=psi_tPlusOne)
                 
                 
    psi_tPlusOneVector <- scoreFunctionUpdate(Z=ZtPlusOne, W=WtPlusOne, SigmaTilde=sigmaTilde, lambdaGrid=lambdaGrid, numKnots=numKnots, gridStep=gridStep,
                                    nFolds=nFolds)
    
    ### Values needed to calculate risk
    unmixedScoreVals_tPlusOne <- psi_tPlusOneVector(ZtPlusOne, derivativeOrder = 0)
    
    unmixedScoreGradient_tPlusOne <- psi_tPlusOneVector(ZtPlusOne, derivativeOrder = 1)
    
    ### The risk at the updated unmixing matrix and score vector. Calculated based on the unmixed score vector.
    risk_tPlusOne <- evaluateRisk(thetas=thetaMatrix, unmixedScoreFunction=psi_tPlusOneVector,
                                  SigmaTilde=sigmaTilde, W=WtPlusOne)
    
    mixedScoreVals_tPlusOne <- t(WtPlusOne) %*% unmixedScoreVals_tPlusOne

    ### The absolute percentage change in the risk
    Change <- (risk_t-risk_tPlusOne)/abs(risk_t)
    
    overMinIterationsU <- iterator > algorithmSpecs["min_iterations"]
    
    closerThanDistU <- Change < algorithmSpecs["min_change"]
    
    if (overMinIterationsU & closerThanDistU){
      
      break
    
    }else{
      
      Wt <- WtPlusOne
      
      Zt <- ZtPlusOne
      
      #psi_t <- psi_tPlusOne
      
      psi_tVector <- psi_tPlusOneVector
      
      unmixedScoreVals_t <- unmixedScoreVals_tPlusOne
      
      mixedScoreVals_t <- mixedScoreVals_tPlusOne
      
      risk_t <- risk_tPlusOne
      
    }
    
    #print(Change)
    
  }
    
  ### Return final U needed for calculating unmixing matrices, and the score vector                    
    
    UandScoreVec <- list(UtPlusOne, psi_tPlusOneVector)
    
    names(UandScoreVec) <- c("U", "score_vector")

  return(UandScoreVec)
  
}
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           


### Update Unmixing Matrix

unmixingMatrixUpdate <- function(thetaStarMatrix, UInit, sigmaTilde,
                                 scoreFunctionVectorInit, minIterationsU=10, maxIterationsU=100,
                                 minChangeU=1e-06, alpha=.2){
    
  Gamma <- cov(t(thetaStarMatrix))
    
  GammaMinusHalf <- solve(sqrtm(Gamma))
    
  thetaTildeMatrix <- GammaMinusHalf %*% thetaStarMatrix
  
  ### Note that A is a symmetric matrix
  A <- GammaMinusHalf %*% sigmaTilde %*% GammaMinusHalf
    
  ### The value of the risk
  rHat <- function(uRisk, ourScoreVector, thetaTildeMatrix, A){ 
      
    UtAUrisk <- uRisk %*% A %*% t(uRisk)
    
    uncoupledThetaTildesUrisk <- uRisk %*% thetaTildeMatrix
      
    riskScoreVals <- ourScoreVector( uncoupledThetaTildesUrisk, derivativeOrder=0)
      
    riskScoreValsPrime <- ourScoreVector( uncoupledThetaTildesUrisk, derivativeOrder=1)
    
    ## This should be trace
    firstPart <- tr(t(riskScoreVals) %*% UtAUrisk %*% riskScoreVals)/dim(riskScoreVals)[2]
      
    secondPart <- mean(2*array(diag(UtAUrisk), dim=c(1, dim(riskScoreValsPrime)[1])) %*% riskScoreValsPrime)
    
    return(firstPart+secondPart)
    
  }
  
  ### The gradient of the risk
  nablaRHat <- function(uRisk, thetaTildeMatrix, scoreFunctionVector, A){      
      
    alternateCrossTermCalculatorForThetaTilde <- function(thetaMat, scoreVectorMat, 
                                                          scoreVectorPrimeMat, UtAU){
        
        finalTermMatrix <- array(NA, dim=rep(dim(scoreVectorMat)[1], 2))
        
        IkByk <- diag(rep(1, dim(scoreVectorMat)[1]))
        
        for (l in 1:dim(thetaTildeMatrix)[1]){
        
            e_l <- array(rep(0, dim(scoreVectorMat)[1]), dim=c(dim(scoreVectorMat)[1], 1) )

            e_l[l,1] <- 1

            neededPrimeVec <- scoreVectorPrimeMat[l,]

            G_l <- scoreVectorMat * sapply(neededPrimeVec, FUN=rep, dim(scoreVectorMat)[1])

            finalTermMatrix[l,] <- (2/dim(scoreVectorMat)[2])*(UtAU[l,] %*% (IkByk-(e_l %*% t(e_l)))
                                                              %*% G_l %*% t(thetaMat)) 
            
            }
        
        return(finalTermMatrix)

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
      
      
    thirdPartU <- (1/dim(thetaTildeTermsMatrix)[2]) * (scoreVals %*% t(scoreVals)-diag(diag(scoreVals %*% t(scoreVals)))) %*% (uRisk %*% t(A))

      
    thirdPartThetaCrossTerms <-  alternateCrossTermCalculatorForThetaTilde(thetaMat=thetaTildeMatrix,
                                 scoreVectorMat=scoreVals, scoreVectorPrimeMat=scoreValsPrime, 
                                  UtAU=UtAU)
      
    nablaValues <- firstPart + secondPart + thirdPartU + thirdPartThetaCrossTerms

    return(nablaValues)
    
  } 
  
  scoreFunctionVector <- scoreFunctionVectorInit
  
  UtMinus1 <- UInit
  
  dist <- 1
  
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
                           
    doesItSatisfy[is.na(doesItSatisfy)] <- TRUE
    
    
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

               
                           
                           
                           
                           
### SCORE FUNCTION WITHOUT SEPARABILITY                           
                           
scoreFunctionUpdate <- function(Z, W, SigmaTilde, lambdaGrid=10^seq(-6, 6, 1), 
                                numKnots=10, gridStep=.01,
                                    nFolds=5){
    
  options(warn=-1)

  basisEvalInLongVec <- function(basisList, zVec){
    
    
    basisEval <- unlist(lapply(1:length(basisList), 
                               FUN=function(x) predict(basisList[[x]], zVec[x])))
                               
    return(array(basisEval, dim=c(length(basisEval), 1)))
    
    }
                               
                               
  evalInBlockDiag <- function(basisList, zVec){
    
    basisEval <- lapply(1:length(basisList), 
                        FUN=function(x) predict(basisList[[x]], zVec[x]))
                               
    return(as.array(bdiag(basisEval)))
    
    }
                        
                        
  evalAsList <- function(basisList, zVec){
    
    basisEval <- lapply(1:length(basisList), 
                        FUN=function(x) predict(basisList[[x]], zVec[x]))
                               
    return(basisEval)
    
    }
                        
                        
  knotMaker <- function(z, numKnots=10){
    
      num_knots_plus_1 <- numKnots+1 ### Give the number of knots, plus 1
  
      theProbs <- seq(1:(num_knots_plus_1-1))/(num_knots_plus_1)
  
      theKnots <- quantile(Z, probs = theProbs) ### Knots based on quantiles
  
  ### Making sure knots are not placed are extreme points in the data
      if (min(theKnots) == min(z)){
    
        theKnots[which(theKnots==min(z))] <- min(z[which(z != min(z))])
    
    
      }
  
      if (max(theKnots) == max(z)){
    
        theKnots[which(theKnots==max(z))] <- min(Z[which(z != max(z))])
    
    
      }
  
  ### Generating the basis to be used to represent the score function
      theKnots <- unique(theKnots)
    
      return(theKnots)
    
    }


    gridMaker <- function(z, gridStep=.001){
    
      zMin <- floor(min(z, na.rm = T))
  
      zMax <- ceiling(max(z, na.rm = T))
  
      zGrid <- seq(zMin, zMax,  gridStep)
    
    }
  
  ### Calculates the risk
  riskCalculator <- function(betaHat, lambda, T1, T2, cW){
    
    riskValue <- t(betaHat) %*% T1 %*% betaHat + 2*t(cW) %*% T2 %*% betaHat
    
    return(riskValue)
    
  }
  
  ### Calculates the basis coefficient vector for the score function
  betaCalculator <- function(lambda, T1, T2, T3, cW){
          
    betaHat <- -1*ginv(T1+lambda*T3) %*% t(T2) %*% cW
      
    
    return(betaHat)
    
  }
  
  ### If bigger than 1, perform cross validation, otherwise estimate from given turning parameter
  CV <- length(lambdaGrid) > 1
    
    
  if (CV & (nFolds==1)){
      
      stop("Must have multiple folds if doing CV.")
      
  }
  
  ### Number of curves
  n <- dim(Z)[2]

  K <- dim(Z)[1]
  
  knots <- alply(Z, 1, knotMaker, numKnots=10)

  grids <- alply(Z, 1, gridMaker)
    
    #### number of knots = df-degree
   #### Specify knots to be at quantiles 
  
  df <- numKnots+1+min(3,n)
                        
  bases<- lapply(1:K,FUN = function(x) bSpline(grids[[x]], degree = min(3,n), df=df,
                                               intercept = T)) ### use knots
                 
                 
  WSigmaTildeTW <- W %*% SigmaTilde %*% t(W)
      
  oneVec <- array(rep(1, dim(W)[1]), dim=c(dim(W)[1], 1))  
      
  cW <- array(apply(W, MARGIN=1, FUN=function(x) t(oneVec) %*% ( SigmaTilde * x %*% t(x) ) %*% oneVec), dim=c(length(oneVec), 1))
      
  n <- dim(Z)[2]
                    
  oneMatrixForKronecker <- array(1, dim=rep(dim(bases[[1]])[2], 2))
      
  WSigmaTildeTWKronecker <-  kronecker(WSigmaTildeTW, oneMatrixForKronecker)
    
  allTheEvals <- alply(Z, 2, basisEvalInLongVec, basisList=bases)
      
  allTheCrossProds <- lapply(allTheEvals, FUN=function(x) crossprod(t(x)))  
      
  basesDerivs <- lapply(bases, FUN=deriv, derivs=1L)
      
  derivEvals <- alply(Z, 2, evalInBlockDiag, basisList=basesDerivs)
                             
  basesSecondDerivs <- lapply(bases, FUN=deriv, derivs=2L)
                             
  T1 <- WSigmaTildeTWKronecker*Reduce('+', allTheCrossProds)/length(allTheCrossProds)
      
  T2 <- Reduce("+", derivEvals)/length(derivEvals) 
    
  T3 <- as.matrix(bdiag(lapply(basesSecondDerivs, FUN=function(x) gridStep*base::crossprod(x) )))
  
  if (CV){ ### Cross-validation for tuning parameter estimation
    
    cvFolds <- createFolds(1:n, k =nFolds)
    
    ### Cross validation done in parallel
    cvErrorsForLambdas <- foreach(lamb=lambdaGrid,
                .packages=c('splines2', "caret", 'MASS', 'plyr', 'Matrix'), .combine = 'c') %dopar% {
      
      risksThisLambda <- c()
      
      cvFolds <- createFolds(1:n, k =nFolds)
      
      for (fold in cvFolds){
      
        trainZ <- Z[, -fold]
        
        valZ <- Z[, fold]
          
        T1Train <- WSigmaTildeTWKronecker*Reduce('+', allTheCrossProds[-fold])/length(allTheCrossProds[-fold])
          
        T1Val <- WSigmaTildeTWKronecker*Reduce('+', allTheCrossProds[fold])/length(allTheCrossProds[fold])
      
        T2Train <- Reduce("+", derivEvals[-fold])/length(derivEvals[-fold]) 
          
        T2Val <- Reduce("+", derivEvals[fold])/length(derivEvals[fold]) 
    
        T3Train <- as.matrix(bdiag(lapply(basesSecondDerivs, FUN=function(x) gridStep*base::crossprod(x) )))
      
        betaHatThisFoldThisLambda <- betaCalculator(lambda=lamb, T1=T1Train, T2=T2Train, T3=T3, cW=cW)
        
        riskThisFoldThisLambda <- riskCalculator(betaHat=betaHatThisFoldThisLambda, 
                                                 lambda=lamb, T1=T1Val, T2=T2Val, cW=cW)
        
        risksThisLambda <- c(risksThisLambda, riskThisFoldThisLambda)
      
      }
     
      mean(risksThisLambda)
       
    }
    
    chosenLambda <- lambdaGrid[which.min(cvErrorsForLambdas)]
    
  }else{
    
    chosenLambda <- lambdaGrid
    
  }
  
  bigBetaHat <- betaCalculator(lambda=chosenLambda, T1=T1, T2=T2, T3=T3, cW=cW)
                                          
  basesDimensionSizes <- unlist(lapply(bases, FUN=function(x) dim(x)[2]))
                                
  betaBreakPoints <- cumsum(basesDimensionSizes)-basesDimensionSizes+1
                                          
  brokenBeta <- Map(function(i,j) bigBetaHat[i:j, , drop=FALSE], 
                    betaBreakPoints, cumsum(diff(c(betaBreakPoints, dim(bigBetaHat)[1]+1))))
                                          
  scoreFunctionCoordinate <- function(basis, betaHat){
      
    finalCoordinate <- function(z, derivativeOrder=0){ ### Note, z is a scalar!
        
        if (derivativeOrder==0){

          basisNow = basis

        }else{

        basisNow <- deriv(basis, derivs = derivativeOrder)
    
    }
    
    basisVec = predict(basisNow, newx = z)
    
    return(basisVec %*% betaHat)
        
        }
      
      return(finalCoordinate)
    
  }
            
  #force(scoreFunctionCoordinate)
                                       
  scoreFunctionList <- lapply(1:length(brokenBeta), FUN = function(k)
      scoreFunctionCoordinate(basis=bases[[k]], betaHat=brokenBeta[[k]])  )      
                              
  scoreListToVector <- function(aScoreList){
    
    ### The output. Its arguments are a numberic vector, and an argument for desired order of derivative
    scoreVector <- function(theta, derivativeOrder=0){
      
      if (!is.matrix(theta)){
        
        theta <- array(theta, dim=c(length(theta), 1))
        
      }
      
      return(t(sapply(1:dim(theta)[1], FUN=function(x) aScoreList[[x]](theta[x,], derivativeOrder=derivativeOrder))))
      
    }
    
    return(scoreVector)
    
    
  }
               
  finalScoreVector <- scoreListToVector(aScoreList = scoreFunctionList)
  
  return(finalScoreVector)
  
}










