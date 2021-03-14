#' @title FCS-Twofold
#' @description An R package for data imputation
#' @field formula Formula
#' @field data A data frame
#' @importFrom methods new
#' @examples
#' data(longitudinal)
#' TwoFCS(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL)
#' TwoFold(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL,...)

TwoFCS <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL) {
  TwoFold(originaldata,smtype,smformula,method,predictorMatrix,m,numit,rjlimit,noisy,errorProneMatrix=errorProneMatrix)
}
TwoFold <- function(originaldata,smtype,smformula,method,predictorMatrix=NULL,m=5,numit=10,rjlimit=1000,noisy=FALSE,errorProneMatrix=NULL,
                        ...) {
  
  dtsamOutcomeDens <- function(inputData) {
    inputDataN <- dim(inputData)[1]
   
    outmodxb <- matrix(outcomeModBeta[1:nTimePoints], nrow=inputDataN, ncol=nTimePoints,byrow=TRUE)
    covXbEffects <-  model.matrix(as.formula(paste("~-1+",strsplit(smformula, "~")[[1]][2],sep="")),
                                  inputData) %*% tail(outcomeModBeta,length(outcomeModBeta)-nTimePoints)
   
    outmodxb <- outmodxb + matrix(covXbEffects, nrow=inputDataN, ncol=nTimePoints)
    
    prob <- Expont(outmodxb)
    logSurvProb <- log(1-prob)
    logSurvProbCumSum <- cbind(rep(0,inputDataN),t(apply(logSurvProb, 1, cumsum)))
    
    lastSurvPlusOne <- 1 + inputData[,timeCol] - inputData[,dCol]
    logSurvProbIndividual <- logSurvProbCumSum[cbind(1:inputDataN,lastSurvPlusOne)]
    exp(logSurvProbIndividual + d*log(prob[cbind(1:inputDataN, originaldata[,timeCol])]))
  }
  
  extraArgs <- list(...)
  
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")
  
  n <- dim(originaldata)[1]
    r <- 1*(is.na(originaldata)==0)
  
  outcomeCol <- which(colnames(originaldata)==as.formula(smformula)[[2]])
  smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]

  partialVars <- which((method=="norm") | (method=="latnorm") | (method=="logreg") | (method=="poisson") | (method=="podds") | (method=="mlogit"))
  if (length(partialVars)==0) stop("You have not specified any valid imputation methods in the method argument.")
  
  #check that methods are given for each partially observed column, and not given for fully observed columns
  for (colnum in 1:ncol(originaldata)) {
    if (method[colnum]!="") {
      #an imputation method has been specified
      if (colnum %in% outcomeCol) {
        stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                   ". Elements of the method argument corresponding to the outcome variable(s) should be empty.",sep=""))
      }
      else {
        if (sum(r[,colnum])==n) {
          stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                     ", but it appears to be fully observed.",sep=""))
        }
      }
    }
    else {
      #no imputation method has been specified
      if (sum(r[,colnum])<n) {
        #some values are missing
        if (((colnum %in% outcomeCol)==FALSE) & (sum(errorProneMatrix[,colnum])==0)) {
          stop(paste("Variable ",colnames(originaldata)[colnum], " does not have an imputation method specified, yet appears to have missing values.",sep=""))
        }
      }
    }
    
  }

  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))
  
  #passive variables
  passiveVars <- which((method!="") & (method!="norm") & (method!="logreg") & (method!="poisson") & (method!="podds") & (method!="mlogit") & (method!="latnorm"))
  
  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))
  
  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }
  
  rjFailCount <- 0
  
  for (imp in 1:m) {
    
    print(paste("Imputation Number:",imp))
    
    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      targetCol <- partialVars[var]
      
      imputations[[imp]][r[,targetCol]==0,targetCol] <- sample(imputations[[imp]][r[,targetCol]==1,targetCol], size=sum(r[,targetCol]==0), replace=TRUE)
      
    }
    
    
    if (sum(r[,outcomeCol])<n) {
      if (imp==1) {
        print("Imputing missing outcomes")
      }
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      
      imputationNeeded <- (1:n)[r[,outcomeCol]==0]
      
      
      ymod <- stats::lm(as.formula(smformula),imputations[[imp]])
      beta <- ymod$coef
      sigmasq <- summary(ymod)$sigma^2
      
      imputations[[imp]][imputationNeeded,outcomeCol] <- 0
      outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
      imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
    }
    
    for (cyclenum in 1:numit) {
      
      if (noisy==TRUE) {
        print(paste("Iteration ",cyclenum))
      }
      #update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[! partialVars %in% targetCol], fullObsVars)
        }
        else {
          predictorCols <- which(predictorMatrix[targetCol,]==1)
          #ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[! predictorCols %in% outcomeCol]
        }
        if ((imp==1) & (cyclenum==1)) {
          
          print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[predictorCols],collapse=',')," plus outcome",collapse=','))
          
        }
        if (length(predictorCols)>0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        }
        else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1",sep=""))
        }
        
        xmoddata <- imputations[[imp]]
        
        if (method[targetCol]=="norm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=xmoddata)
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          
          xfitted <- model.matrix(xmod) %*% newbeta
          
        }  
        if (noisy==TRUE) {
          print(summary(xmod))
        }
        
        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          varcov <- vcov(ymod)
          outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
          covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
          outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        
        if ((imp==1) & (cyclenum==1) & (var==1)) {
          
          smCoefIter <- array(0, dim=c(m, length(outcomeModBeta), numit))
          
        }
        
        if (var==length(partialVars)) {
          #then we have reached end of a cycle
          
          smCoefIter[imp,,cyclenum] <- outcomeModBeta
          
        }
        
        #impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[,targetCol]==0]
        
        if ((method[targetCol]=="logreg") | (method[targetCol]=="podds") | (method[targetCol]=="mlogit")) {
          #directly sample
          if (method[targetCol]=="logreg") {
            numberOutcomes <- 2
            fittedMean <- cbind(1-xfitted, xfitted)
          }
          else {
            numberOutcomes <- nlevels(imputations[[imp]][,targetCol])
            fittedMean <- xfitted
          }
          
          outcomeDensCovDens = array(dim=c(length(imputationNeeded),numberOutcomes),0)
          
          for (xMisVal in 1:numberOutcomes) {
            if (method[targetCol]=="logreg") {
              if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
                valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
              }
              else {
                valToImpute <- xMisVal-1
              }
            }
            else {
              valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
            }
            imputations[[imp]][imputationNeeded,targetCol] <- valToImpute
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              outcomeDens <- dnorm(deviation, mean=0, sd=outcomeModResVar^0.5)
            } 
            
            outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
          }
          directImpProbs = outcomeDensCovDens / rowSums(outcomeDensCovDens)
          
          if (method[targetCol]=="logreg") {
            directImpProbs = directImpProbs[,2]
            if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
              imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[1]
              imputations[[imp]][imputationNeeded,targetCol][rbinom(length(imputationNeeded),1,directImpProbs)==1] <- levels(imputations[[imp]][,targetCol])[2]
            }
            else {
              imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),1,directImpProbs)
            }
          }
          else {
            imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[apply(directImpProbs, 1, catdraw)]
          }
          
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        else {
          #use rejection sampling
          #first draw for all subjects who need imputing, using a small number of attempts
          firstTryLimit <- 25
          j <- 1
          
          while ((length(imputationNeeded)>0) & (j<firstTryLimit)) {
            #sample from covariate model
            if ((method[targetCol]=="norm") | (method[targetCol]=="latnorm")) {
              imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
            }
            else if (method[targetCol]=="poisson") {
              imputations[[imp]][imputationNeeded,targetCol] <- rpois(length(imputationNeeded),xfitted[imputationNeeded])
            }
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            #accept/reject
            uDraw <- runif(length(imputationNeeded))
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(length(imputationNeeded),1))))
            }
            
            imputationNeeded <- imputationNeeded[reject==1]
            
            j <- j+1
          }
          
          #now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {
            
            tempData <- imputations[[imp]][i,]
            tempData <- tempData[rep(1,rjlimit),]
            if (method[targetCol]=="norm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq^0.5)
            }
            else if (method[targetCol]=="logreg") {
              tempData[,targetCol] <- rbinom(rjlimit,size=1,xfitted[i])
            }
            else if (method[targetCol]=="poisson") {
              tempData[,targetCol] <- rpois(rjlimit,xfitted[i])
            }
            else if (method[targetCol]=="latnorm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq[i]^0.5)
            }
            
            #passively impute
            tempData <- updatePassiveVars(tempData, method, passiveVars)
            
            #accept reject
            uDraw <- runif(rjlimit)
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              deviation <- tempData[,outcomeCol] - outmodxb
              reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(rjlimit,1))))
            }
            
            if (sum(reject)<rjlimit) {
              imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
            } else {
              rjFailCount <- rjFailCount + 1
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
      }
      
      #imputations of missing outcomes, if present (using proper imputation), for regression and logistic
      #substantive models
      if ((smtype=="lm") ) {
        if (sum(r[,outcomeCol])<n) {
          imputationNeeded <- (1:n)[r[,outcomeCol]==0]
          #estimate parameters of substantive model using those with outcomes observed
          if (smtype=="lm") {
            ymod <- lm(as.formula(smformula),imputations[[imp]][r[,outcomeCol]==1,])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
            covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
            outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
            imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
          }
          
        }
      }
    }
    
  }
  
  if (rjFailCount>0) {
    warning(paste("Rejection sampling failed ",rjFailCount," times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.",sep=""))
  }
  
  # Added smformula and smtype to metadata, and make "smcfcs class"
  res <- list(
    impDatasets = imputations,
    smCoefIter = smCoefIter,
    smInfo = list("smtype" = smtype, "smformula" = smformula)
  )
  class(res) <- "smcfcs"
  
  return(res)
}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

Expont <- function(x) {
  exp(x)/(1+exp(x))
} 

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1,size=1,prob=prob)==1]
}

#########################################################
#########################################################


