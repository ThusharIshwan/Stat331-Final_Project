setwd("C:/Users/Thushar/Desktop/Stat331/Project")
# Outlier Take out Simulation

#What data are we using:
load("pollution.Rdata")

y = pollution$e3_bw
X = pollution[!names(pollution) %in% c("e3_bw")]


# Splitting the columns into four domains
pollution_chemicals <- pollution[c(18:27,31:70)] # chemical domain
pollution_outdoors <- pollution[c(2:5,28:30,71:73)] # outdoor exposures domain
pollution_lifestyles <- pollution[c(6:17)] # lifestyles domain
pollution_others <- pollution[c(74:80)] # covariates domain


# A function for finding the DFFits outliers at a specific tolerance. Ploting the values is an included option.
Get_DFFITS_Outliers <- function(M, tol = 2, to_plot = FALSE, ylab = "DFFITS")
{
  D = dffits(M)
  n = nobs(M)
  p <- length(M$coef)-1
  Out_Indices = which(abs(D)>tol*sqrt((p+1)/n))
  
  if (to_plot)
  {
    plot (D, ylab = ylab)
    abline(h=tol*sqrt((p+1)/n),lty=2)
    abline(h=-tol*sqrt((p+1)/n),lty=2)
    
    points(D[Out_Indices]~Out_Indices,col="red",pch=19)
    text(y=D[Out_Indices],x=Out_Indices, labels=Out_Indices, pos=2)
    
  }
  
  return(as.vector(Out_Indices))
}

# One instance of the simulation
Simulation_Instance <- function(y_train, X_train, y_test, X_test, cols, plot_outliers = FALSE)
{
  
  M = lm(y_train~.,data = X_train[cols]) # Create the model using the entire training set
  M_Outliers = c()
  
  M_Outliers = Get_DFFITS_Outliers(M, to_plot = plot_outliers) # Identify the outliers 
  
  
  X_train_rm_out = X_train[-M_Outliers,cols]
  y_train_rm_out = y_train[-M_Outliers]
  
  M_Outliers_Removed = lm(y_train_rm_out~.,data =X_train_rm_out) # Remove the outliers from the training set and remake the model
  
  mspe_1 = MSPE(y_test, X_test, M)
  mspe_2 = MSPE(y_test, X_test, M_Outliers_Removed) #Compare the MSPEs of the 2 models
  
  return(list(length(cols),length(M_Outliers),mspe_1, mspe_2, mspe_1-mspe_2)) # Return information about the instance
  
}

# Full simulation:
Run_Simulation <- function(ITERS, y, X, CV = 1, tts = 0.6, prop0 = 1, prop1 = 1)
{
  COLS = length(X) # number of columns
  table = data.frame(features = rep(NA,ITERS*CV),num_outliers = rep(NA,ITERS*CV),mspe_original = rep(NA,ITERS*CV),
                     mspe_no_outliers = rep(NA,ITERS*CV),mspe_difference = rep(NA,ITERS*CV))
  
  for (cv in 1:CV)
  {
    sampler = 1:(length(y) * tts) + (((length(y)*(1-tts) * (cv - 1)) / (CV - 1))) # different training splits used for cross-validation
    
    y_train = y[sampler]
    X_train = X[sampler,]
    y_test = y[-sampler]
    X_test = X[-sampler,]
    
    # Creates the bool matrix of chosen columns
    RandomFeatureSelectionList <- list()
    for (i in 1:ITERS)
    {
      RandomFeatureSelectionList[[i]] = sample(c(rep(0,prop0),rep(1,prop1)), replace=TRUE, size=COLS)
    }
    # Converts The booleans of the matrix into a list of columns
    
    for (i in 1:length(RandomFeatureSelectionList))
    {
      col_list = c()
      for(j in 1:COLS)
      {
        if(RandomFeatureSelectionList[[i]][j] == 1)
        {
          col_list = append(col_list,j)
        }
      }
      sim_inst = try(Simulation_Instance(y_train, X_train, y_test, X_test, col_list), silent = TRUE)
      if (length(sim_inst) == 5)
      {
        for (j in 1:5)
        {
          
          table[((cv-1)*ITERS) + i,j] = sim_inst[[j]]
          
        }
      }
    }

  }
  
  return(table)
}

# Draws a histogram and the ascosiated normal distribution
Hist_Norm <- function(d, breaks = 20, x_lab = "x_lab", y_lab = "y_lab", main = "main"){
  full = d[!is.na(d)]
  print(min(full))
  h <- hist(full, breaks = breaks, x_lab = "x_lab", y_lab = "y_lab", main = "main")
  xfit <- seq(min(full), max(full), length = 40) 
  yfit <- dnorm(xfit, mean = mean(full), sd = sd(full))
  yfit <- yfit * diff(h$mids[1:2]) * length(full)
  lines(xfit, yfit, col = "black", lwd = 2)
  
}


# Set the seed for ease of replication:
set.seed(20776408)

# Models with equal chance of including or excluding columns, across the different domains
domain1 = Run_Simulation(1000,y,pollution_chemicals, CV = 5)
domain2 = Run_Simulation(1000,y,pollution_outdoors, CV = 5)
domain3 = Run_Simulation(1000,y,pollution_lifestyles, CV = 5)
domain4 = Run_Simulation(1000,y,pollution_others, CV = 5)
domain12 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors), CV = 5)
domain13 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles), CV = 5)
domain14 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_others), CV = 5)
domain23 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles), CV = 5)
domain24 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_others), CV = 5)
domain34 = Run_Simulation(1000,y,cbind(pollution_lifestyles, pollution_others), CV = 5)
domain123 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles), CV = 5)
domain124 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_others), CV = 5)
domain134 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles, pollution_others), CV = 5)
domain234 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5)
domain1234 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5)

# Models with a higher chance of excluding columns, across the different domains
low_col_domain1 = Run_Simulation(1000,y,pollution_chemicals, CV = 5, prop0 = 2)
low_col_domain2 = Run_Simulation(1000,y,pollution_outdoors, CV = 5, prop0 = 2)
low_col_domain3 = Run_Simulation(1000,y,pollution_lifestyles, CV = 5, prop0 = 2)
low_col_domain4 = Run_Simulation(1000,y,pollution_others, CV = 5, prop0 = 2)
low_col_domain12 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors), CV = 5, prop0 = 2)
low_col_domain13 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles), CV = 5, prop0 = 2)
low_col_domain14 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_others), CV = 5, prop0 = 2)
low_col_domain23 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles), CV = 5, prop0 = 2)
low_col_domain24 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_others), CV = 5, prop0 = 2)
low_col_domain34 = Run_Simulation(1000,y,cbind(pollution_lifestyles, pollution_others), CV = 5, prop0 = 2)
low_col_domain123 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles), CV = 5, prop0 = 2)
low_col_domain124 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_others), CV = 5, prop0 = 2)
low_col_domain134 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles, pollution_others), CV = 5, prop0 = 2)
low_col_domain234 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5, prop0 = 2)
low_col_domain1234 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5, prop0 = 2)

# Models with a higher chance of including columns, across the different domains
high_col_domain1 = Run_Simulation(1000,y,pollution_chemicals, CV = 5, prop1 = 2)
high_col_domain2 = Run_Simulation(1000,y,pollution_outdoors, CV = 5, prop1 = 2)
high_col_domain3 = Run_Simulation(1000,y,pollution_lifestyles, CV = 5, prop1 = 2)
high_col_domain4 = Run_Simulation(1000,y,pollution_others, CV = 5, prop1 = 2)
high_col_domain12 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors), CV = 5, prop1 = 2)
high_col_domain13 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles), CV = 5, prop1 = 2)
high_col_domain14 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_others), CV = 5, prop1 = 2)
high_col_domain23 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles), CV = 5, prop1 = 2)
high_col_domain24 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_others), CV = 5, prop1 = 2)
high_col_domain34 = Run_Simulation(1000,y,cbind(pollution_lifestyles, pollution_others), CV = 5, prop1 = 2)
high_col_domain123 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles), CV = 5, prop1 = 2)
high_col_domain124 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_others), CV = 5, prop1 = 2)
high_col_domain134 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_lifestyles, pollution_others), CV = 5, prop1 = 2)
high_col_domain234 = Run_Simulation(1000,y,cbind(pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5, prop1 = 2)
high_col_domain1234 = Run_Simulation(1000,y,cbind(pollution_chemicals, pollution_outdoors, pollution_lifestyles, pollution_others), CV = 5, prop1 = 2)


#combining the observations across the column inclusion probabilities
tot_domain1 = unique(rbind(low_col_domain1, domain1, high_col_domain1))
tot_domain2 = unique(rbind(low_col_domain2, domain2, high_col_domain2))
tot_domain3 = unique(rbind(low_col_domain3, domain3, high_col_domain3))
tot_domain4 = unique(rbind(low_col_domain4, domain4, high_col_domain4))
tot_domain12 = unique(rbind(low_col_domain12, domain12, high_col_domain12))
tot_domain13 = unique(rbind(low_col_domain13, domain13, high_col_domain13))
tot_domain14 = unique(rbind(low_col_domain14, domain14, high_col_domain14))
tot_domain23 = unique(rbind(low_col_domain23, domain23, high_col_domain23))
tot_domain24 = unique(rbind(low_col_domain24, domain24, high_col_domain24))
tot_domain34 = unique(rbind(low_col_domain34, domain34, high_col_domain34))
tot_domain123 = unique(rbind(low_col_domain123, domain123, high_col_domain123))
tot_domain124 = unique(rbind(low_col_domain124, domain124, high_col_domain124))
tot_domain134 = unique(rbind(low_col_domain134, domain134, high_col_domain134))
tot_domain234 = unique(rbind(low_col_domain234, domain234, high_col_domain234))
tot_domain1234 = unique(rbind(low_col_domain1234, domain1234, high_col_domain1234))



#Calcluating mean and standard deviation of the MSPE differences observed above, and plotting the respective histograms.
Hist_Norm(tot_domain1234$mspe_difference, xlab = "Original Model MSPE minus Outlier Removed Model MSPE", main = "MSPE Difference in models across all domains")
mean(tot_domain1234$mspe_difference, na.rm = TRUE)
sd(tot_domain1234$mspe_difference, na.rm = TRUE)

Hist_Norm(tot_domain1$mspe_difference, xlab = "Original Model MSPE minus Outlier Removed Model MSPE", main = "MSPE Difference in models within the Chemicals Domain")
mean(tot_domain1$mspe_difference, na.rm = TRUE)
sd(tot_domain1$mspe_difference, na.rm = TRUE)
Hist_Norm(tot_domain2$mspe_difference, xlab = "Original Model MSPE minus Outlier Removed Model MSPE", main = "MSPE Difference in models Within the Outdoors Domain")
mean(tot_domain2$mspe_difference, na.rm = TRUE)
sd(tot_domain2$mspe_difference, na.rm = TRUE)
Hist_Norm(tot_domain3$mspe_difference, xlab = "Original Model MSPE minus Outlier Removed Model MSPE", main = "MSPE Difference in models within the Lifestyle Domain")
mean(tot_domain3$mspe_difference, na.rm = TRUE)
sd(tot_domain3$mspe_difference, na.rm = TRUE)
Hist_Norm(tot_domain4$mspe_difference, xlab = "Original Model MSPE minus Outlier Removed Model MSPE", main = "MSPE Difference in models Within the Misc. Domain")
mean(tot_domain4$mspe_difference, na.rm = TRUE)
sd(tot_domain4$mspe_difference, na.rm = TRUE)
