setwd("~/Documents/R/STAT331")
load("/Users/harry/Documents/R/STAT331/pollution.Rdata")
pollution$e3_alcpreg_yn_None <- as.factor(pollution$e3_alcpreg_yn_None)
pollution$h_folic_t1_None <- as.factor(pollution$h_folic_t1_None)
pollution$e3_yearbir_None <- as.factor(pollution$e3_yearbir_None)
pollution$h_edumc_None <- as.factor(pollution$h_edumc_None)

e3_bw <- pollution$e3_bw # response
pollution_chemicals <- select(pollution, c(18:27,31:70)) # chemical domain
pollution_outdoors <- select(pollution, c(2:5,28:30,71:73)) # outdoor exposures domain
pollution_lifestyles <- select(pollution, c(6:17)) # lifestyles domain
pollution_others <- select(pollution, c(74:80)) # covariates domain

data_set <- data.frame(e3_bw,pollution_chemicals)

# splitting training set and test set
training_set <- data_set[c(1:600),]
test_set <- data_set[-c(1:600),]

# Forward Selection
M_base <- lm(e3_bw ~ 1, data=training_set)
M_full <- lm(e3_bw ~ (.)^2, data=training_set)
M_forward <- step(object=M_base, scope=list(lower=M_base,upper=M_full),
                  direction="forward", trace=0)

# DFFITS
n <- nobs(M_forward)
p <- length(M_forward$coef)-1
dffits_fwd <- dffits(M_forward) 

## plot DFFITS
plot(dffits_fwd,ylab="DFFITS") 
abline(h=2*sqrt((p+1)/n),lty=2)  ## add thresholds
abline(h=-2*sqrt((p+1)/n),lty=2)
## highlight influential points
dff_ind_fwd <- which(abs(dffits_fwd)>2*sqrt((p+1)/n)) # 33 outliers
points(dffits_fwd[dff_ind_fwd]~dff_ind_fwd,col="red",pch=19) ## add red points
text(y=dffits_fwd[dff_ind_fwd],x=dff_ind_fwd, labels=dff_ind_fwd, pos=2) ## label high influence points

# create training data removing outliers found using DFFITS above
training_set_no_outlier_fwd <- training_set[-dff_ind_fwd,]

# fit model on this training set
M_base_no_outlier <- lm(e3_bw ~ 1, data=training_set_no_outlier_fwd)
M_full_no_outlier <- lm(e3_bw ~ (.)^2, data=training_set_no_outlier_fwd)
M_forward_no_outlier <- step(object=M_base_no_outlier,
                             scope=list(lower=M_base_no_outlier,
                                        upper=M_full_no_outlier),
                             direction="forward", trace=0)

# we use MSPE to assess the prediction accuracy of of the models
MSPE <- function(yi_new, yi_new_hat){
  mean((yi_new - yi_new_hat)^2)
}
yi_new <- test_set$e3_bw

# MSPE on model 1
M1.pred <- predict(M_forward, newdata = test_set[,-1])
MSPE_M1 <- MSPE(yi_new,M1.pred)
MSPE_M1

# MSPE on model 2
M2.pred <- predict(M_forward_no_outlier, newdata = test_set[,-1])
MSPE_M2 <- MSPE(yi_new,M2.pred)
MSPE_M2