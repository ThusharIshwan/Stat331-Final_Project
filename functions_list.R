#Functions for the project



# Methods of Finding Outliers:

Get_Cook_Outliers <- function(M, tol = 0.1, to_plot = FALSE, ylab = "Cook's Distance")
{
  C = cooks.distance(M)
  n = nobs(M)
  p <- length(M$coef)-1
  Out_Indices = which(pf(C,p+1,n-p-1,lower.tail = TRUE) > tol)
  
  if(to_plot)
  {
    plot(M, which = 4)
    plot(C, ylab = ylab)
    points(C[Out_Indices]~Out_Indices,col="red",pch=19)
    text(y=C[Out_Indices],x=Out_Indices, labels=Out_Indices, pos=4)
  }
  
  return(as.vector(Out_Indices))
}


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

# Evaluation Criteria:

MSPE <- function(y, X, M) {
  mean((y - predict(M,X))^2)
}

