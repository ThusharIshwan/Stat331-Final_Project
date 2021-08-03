setwd("C:/Users/Thushar/Desktop/Stat331/Project")
# Outlier Take out Simulation

#What data are we using:
load("pollution.Rdata")

y = pollution$e3_bw
X = pollution[!names(pollution) %in% c("e3_bw")]

pollution_chemicals <- pollution[c(18:27,31:70)] # chemical domain
pollution_outdoors <- pollution[c(2:5,28:30,71:73)] # outdoor exposures domain
pollution_lifestyles <- pollution[c(6:17)] # lifestyles domain
pollution_others <- pollution[c(74:80)] # covariates domain


# One instance of the simulation
Simulation_Instance <- function(y_train, X_train, y_test, X_test, cols, plot_outliers = FALSE)
{
  
  M = lm(y_train~.,data = X_train[cols])
  M_Outliers = c()
  
  M_Outliers = Get_DFFITS_Outliers(M, to_plot = plot_outliers)
  
  
  X_train_rm_out = X_train[-M_Outliers,cols]
  y_train_rm_out = y_train[-M_Outliers]
  
  M_Outliers_Removed = lm(y_train_rm_out~.,data =X_train_rm_out)
  
  mspe_1 = MSPE(y_test, X_test, M)
  mspe_2 = MSPE(y_test, X_test, M_Outliers_Removed)
  
  return(list(length(cols),length(M_Outliers),mspe_1, mspe_2, mspe_1-mspe_2))
  
}

# Full simulation:

Run_Simulation <- function(ITERS, y, X, CV = 1, tts = 0.6, prop0 = 1, prop1 = 1)
{
  COLS = length(X)
  table = data.frame(features = rep(NA,ITERS*CV),num_outliers = rep(NA,ITERS*CV),mspe_original = rep(NA,ITERS*CV),
                     mspe_no_outliers = rep(NA,ITERS*CV),mspe_difference = rep(NA,ITERS*CV))
  
  for (cv in 1:CV)
  {
    sampler = 1:(length(y) * tts) + (((length(y)*(1-tts) * (cv - 1)) / (CV - 1)))
    
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

# Set the seed for ease of replication:
set.seed(20776408)
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

for (r in c(tot_domain1$mspe_difference, 
            tot_domain2$mspe_difference, 
            tot_domain3$mspe_difference, 
            tot_domain4$mspe_difference, 
            tot_domain12$mspe_difference, 
            tot_domain13$mspe_difference, 
            tot_domain14$mspe_difference, 
            tot_domain23$mspe_difference, 
            tot_domain24$mspe_difference, 
            tot_domain34$mspe_difference, 
            tot_domain123$mspe_difference, 
            tot_domain124$mspe_difference, 
            tot_domain134$mspe_difference, 
            tot_domain234$mspe_difference, 
            tot_domain1234$mspe_difference)){
  hist(r)
  
}

