setwd("~/Documents/R/STAT331")
load("/Users/harry/Documents/R/STAT331/pollution.Rdata")
pollution$e3_alcpreg_yn_None <- as.factor(pollution$e3_alcpreg_yn_None)
pollution$h_folic_t1_None <- as.factor(pollution$h_folic_t1_None)
pollution$e3_yearbir_None <- as.factor(pollution$e3_yearbir_None)
pollution$h_edumc_None <- as.factor(pollution$h_edumc_None)

M <- lm(e3_bw ~., data=pollution) # full model no interactions
p <- length(coef(M))-1 # num of covariates
n <- nobs(M) # full observations
X <- model.matrix(M)[,-1] # covariates
y <- pollution$e3_bw # response

# take first 600 to be training
X_train <- X[1:600,] 
X_test <- X[-c(1:600),]
y_train <- y[1:600]
y_test <- y[-c(1:600)]

# lasso model
M_lasso <- cv.glmnet(x=X_train,y=y_train,alpha = 1)

#lambda.best <- M.lasso$lambda.min
#M.lasso <- glmnet(x=model.matrix(M)[,-1],y=data.train.full$e3_bw,alpha = 1,lambda=lambda.best)
#coef(M_lasso, s = "lambda.min")

# prediction in test sset of lasso model with best lambda
pred_lasso <- predict(M_lasso, newx=X_test,  s="lambda.min")

# MSPE
MSPE_lasso <- mean((pred_lasso-y_test)^2)