#' @title Linear Regression
#' @description You can have Reference Class containing some calculations by giving formula and data.
#' @field formula Formula
#' @field data A data frame
#' @importFrom methods new
#' @export ridgereg
#' @exportClass ridgereg
ridgereg <- setRefClass(
  Class = "ridgereg",
  fields = list(
    formula="formula", 
    Fits="numeric",
    Coef="numeric",
    data="data.frame",
    beta_zero = "numeric",
    Call="character",
    beta_hat_reg = "matrix",
    lambda="numeric",mat="matrix"
  ),
  methods = list(
    initialize = function(formula = formula, data = data,lambda=0){
      
      x<-model.matrix(formula,data)
      x <- x[,-1] # Removes (Intercept)
      x[,1:ncol(x)]<-apply(x[,1:ncol(x)],2,function(a) (a-mean(a))/sd(a))
      # x <- scale(x)
      
      qr_mat <- qr(x)
      Q <- qr.Q(qr_mat)
      R <- qr.R(qr_mat)
      
      y<-all.vars(formula)[1]
      y<-as.matrix(data[,names(data)==y])
      temp_data <- data.frame(cbind(x[,2:ncol(x)],y))
      names(temp_data)[ncol(temp_data)]<-all.vars(formula)[1]
      data <<- temp_data
      
      I_lambda<-matrix(nrow=ncol(x),ncol=ncol(x),data = 0)
      
      # Ridge regression coefficients
      beta_hat_reg <<- 
        solve(t(R) %*% R + 0.5 * diag(dim(t(R) %*% R)[1])) %*% t(x) %*% y
      
      b_hat<-(solve((t(x)%*%x)+lambda))%*%(t(x)%*%y)
      y_fits<-x%*%b_hat
      
      coef<-as.numeric(b_hat)
      names(coef)<-rownames(b_hat)
      
      y_fits<-as.numeric(y_fits)
      names(y_fits)<-rownames(y_fits)
      
      beta_zero <<- mean(y)
      Fits <<- y_fits
      Coef <<- coef
      formula<<-formula
      data<<-data
    },
    print = function(){
      "This function prints regression coefficients by using given formula and data in initialization."
      cat("Call:",sep="\n")
      cat(paste("ridgereg(","formula = ",formula[2]," ",formula[1]," ",formula[3],", ","data = ",deparse(substitute(data)),")",sep=""), sep="\n")
      cat(sep="\n")
      cat("Coef:")
      cat(sep="\n")
      
      beta<-Coef
      namn<-names(beta)
      names(beta)<-NULL
      beta<-round(beta,4)
      
      for(i in 2:length(beta)){
        beta[i]<-format(beta[i], width=max(nchar(beta[i]),nchar(namn[i])),justify = "right")
      }
      
      beta[1]<-format(beta[1], width=max(nchar(beta[1]),nchar(namn[1]),nchar("Coef")),justify = "right")
      namn[1]<-format(namn[1], width=max(nchar(beta[1]),nchar(namn[1]),nchar("Coef")),justify = "right")
      
      beta[1]<-paste(beta[1],"  ",sep="")
      namn[1]<-paste(namn[1],"  ",sep="")
      
      beta[2]<-paste(beta[2]," ",sep="")
      namn[2]<-paste(namn[2]," ",sep="")
      
      cat(" ")
      cat(namn)
      cat(" ")
      cat(sep="\n")
      cat(beta)
      
    },
    predict = function(newdata = NULL) {
      "Returns prediction"
      if(is.null(newdata)){
        result <- (Fitted_values = round(Fits, 2))
      } else{
        newdata <- data.frame(newdata)
        X <- as.matrix(scale(newdata))
        beta_final <-matrix(Coef, nrow=length(Coef))
        pred <- (X %*% beta_final) + beta_zero
        result <- pred[,1]
      }
      return(result)
      
    },
    coef = function() {
      coef_v <- as.vector(beta_hat)
      names(coef_v) <- c(row.names(beta_hat))
      return(coef_v)
    }
    
  )
)

# library (MASS)
# lmR <- lm.ridge(Petal.Length ~ Species, iris)
# round(as.vector(lmR$coef),2)
# 
# ridgereg <-  ridgereg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data=iris)
# round(as.vector(ridgereg$Coef[2:3]),2)
# ridgereg$print()
# ridgereg$predict()
# ridgereg$coef()


# library(caret)
# library(mlbench)
# 
# data("BostonHousing")
# boston <- BostonHousing[,-4]
# colnames(boston[-10])
# 
# inTrain <- createDataPartition(y = boston$tax,
#                                p = .8,
#                                list = FALSE)
# training <- boston[ inTrain,]
# testing <- boston[-inTrain,]
# 
# Ridgereg_model <- list(
#   type = c("Regression"),
#   library = "Lab7",
#   loop = NULL,
#   prob = NULL)
# 
# Ridgereg_model$parameters <- data.frame(parameter = "lambda",
#                                         class = "numeric",
#                                         label = "lambda")
# 
# Ridgereg_model$grid <- function (x, y, len = NULL, search = "grid"){
#   data.frame(lambda = seq(0,2,by=0.25))
# }
# 
# Ridgereg_model$fit <- function (x, y, wts, param, lev, last, classProbs, ...){
#   dat <- if (is.data.frame(x))
#     x
#   else
#     as.data.frame(x)
#   dat$.outcome <- y
#   
#   output <-
#     ridgereg$new(.outcome ~ ., data = dat, lambda = param$lambda, ...)
#   
#   output
# }
# 
# 
# Ridgereg_model$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
#   if (!is.data.frame(newdata))
#     newdata <- as.data.frame(newdata)
#   modelFit$predict(newdata)
# }
# set.seed(1991)
# ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
# ridgereg_fit <- caret::train(tax ~ ., data=training, method=Ridgereg_model, trControl=ctrl)
# ridgereg_fit

# ridgereg <-  ridgereg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data=iris)
# ridgereg$predict(train)
# lm_fit <- train(tax ~ .,
#                 data = training,
#                 method = "lm",
#                 trControl = ctrl,
#                 preProc = "scale")
# 
# Test_lm <- predict(lm_fit, testing)
# postResample(Test_lm,testing$tax)



# 
# ridge <- train(lstat ~., data = train,
#                method='leapForward')
# 
# summary(ridge)
# 
# 
# ldaModelInfo <- getModelInfo(model = "leapForward", regex = FALSE)[[1]]
# 
# customModel <- list(label = "Custom Ridgereg implementation", type = "Regression", fit = )
