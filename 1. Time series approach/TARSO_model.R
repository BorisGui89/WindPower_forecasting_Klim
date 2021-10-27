# Clear environment
rm(list = ls())

# Load data
df <- read.csv('C:\\Users\\33637\\Documents\\12. Advanced time series\\11. Exercises\\4. Computer exercise 4\\1. Analyse\\1. Time series approach\\cex4WindDataInterpolated.csv', header = TRUE)

# Select column for calculation
p = as.vector(t(unlist(df[,3, drop=FALSE])))
wd = as.vector(t(unlist(df[,5, drop=FALSE])))

# Number of samples
n = length(p) #changed to 6700 for a series without NaN values

#### Parameters estimation ####

# Loss function creation
TARSO <- function(theta)
{
  # Calculate the estimations based on known data
  p_hat <- numeric(n)
  p_hat[1] <- p[1]
  for (t in 2:n) {
    # Regime selection
    if (!is.na(p[t-1]) & !is.na(wd[t])){
      if (wd[t]>300 | wd[t]<=60){
        p_hat[t] = theta[1] + theta[2]*p[t-1]
      }
      if (wd[t]>60 & wd[t]<=180){
        p_hat[t] = theta[3] + theta[4]*p[t-1]
      }
      if (wd[t]>180 & wd[t]<=300){
        p_hat[t] = theta[5] + theta[6]*p[t-1]
      }
    }
    else{
      p_hat[t] <- NA
    }
  }
  ## Calculate the objective function value
  
  QN <- sum((p-p_hat)^2,na.rm=TRUE)
  
  return(QN)
}

# Test the function
theta <- c(0,1,0,1,0,1)
TARSO(theta)

# Find the minimum of the function
opt = optim(c(0,0,0,0,0,0),TARSO)

RMSE = sqrt(opt$value/n)
CI95 = RMSE*qnorm(0.975)
print(CI95)

params = opt$par
print(params)
