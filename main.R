#############################################################
####### Code for logit model
#############################################################
library(MASS)
source('functions.R')
# 1.generate data------------
## parameters values
lam_real = 1 # Box(X,lam_real) is normally distributed
mean_box_x = 5 # true mean of Box(X,lam_real)
s2_box_x = 1 # true variance of Box(X,lam_real)
s2_u = 0.8 # true variance of U
beta_real = c(-0.42,log(1.5)) # c(beta_0,beta1)
## sample size
n = 1000 
r = 2
## generate dataset
set.seed(123)
X_int = Inv_Box(rnorm(n, mean = mean_box_x, sd = sqrt(s2_box_x)), lam_real) # X has no replicates
U_int = matrix(rnorm(n*r, 0, sqrt(s2_u)), n, r)
W_int = matrix(0, n, r)
for(j in 1:r){W_int[,j] = Inv_Box(Box(X_int,lam_real) + U_int[,j], lam_real)} 
w = W_int[,1] 
## generate response y
pr = H_fn(beta_real[1] + beta_real[2]*X_int)
y = c()
for(i in 1:n){y[i] = rbinom(1,size = 1,prob = pr[i])}
## define cut points
J = 5 # categorize W into five categories based on quintiles
C = rep(0, J-1)
for(j in 1:(J-1)){C[j] = quantile(w,j/J)} 

# 2.estimate theta_J-theta1 using our simfex method------------
## input
w = w # observations with measurement error, a n-dim vector
y = y # response, a n-dim vector
B.boot = 50 # Bootstrap number for variance estimation
## run our simfex method
fn_simfex(w, y, B.boot)
## output
## the point estimator, variance estimator and p-value for theta_J-theta_1
