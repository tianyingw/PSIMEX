#############################################################
###### functions for logit model 
#############################################################
## Logistic function
H_fn <- function(x){
  return(1/(1+exp(-x)))
} 
# Logistic inverse transformation
H_fn_inv <- function(x){
  return(log(x/(1-x)))
}
## Box-Cox transformation：require X>0
Box <- function(x,lam){
  if(lam==0){return(log(x))}else{
    return((x^lam-1)/lam)}}
Inv_Box <- function(y,lam){
  if(lam==0){return(exp(y))}else{
    return((lam*y+1)^(1/lam))}} 
## estimate misclassification matrix and p using the BoxCox method
fn_box <- function(mean_box_x,s2_box_x,s2_u,lam,a){
  fn_x <-  function(x){
    out = dnorm(Box(x,lam),mean=mean_box_x,sd=sqrt(s2_box_x))*x^(lam-1) 
    return(out)}  # f(x)
  fn_wx <- function(w,x){
    out = dnorm(Box(w,lam),mean=Box(x,lam),sd=sqrt(s2_u))*w^(lam-1)
    return(out)}  # f(w|x)
  # Double integrals(w in [a1,a2], x in [b1,b2])
  myfun = function(w,x) { 
    out = fn_wx(w,x)*fn_x(x)
    return(out)}
  fn_d_int <- function(a1,a2,b1,b2){
    out = integrate(function(x) { 
      sapply(x, function(x) {
        integrate(function(w) myfun(w,x), a1, a2)$value
      })
    }, b1, b2)$value
    return(out)}
  # Single integral
  fn_s_int <- function(b1,b2){
    out = integrate(fn_x, b1, b2)$value
    return(out)}
  # estimate transformation matrix A_box: using f(w|x) f(x)
  A_box = matrix(NA,J,J) #w in [a1,a2], x in [b1,b2]
  p_hat = c()
  for(i in 1:J){
    p_hat[i] = fn_s_int(a[i],a[i+1])
    for(j in 1:J){
      A_box[i,j] = fn_d_int(a[j],a[j+1],a[i],a[i+1])/p_hat[i]}}
  A_box = fn_norm(A_box) 
  out = list()
  out$A_box = A_box
  out$p_hat = p_hat
  return(out)
}
## categorized covariate function of x: x is a real value.
M <- function(x){
  Mx = vector()
  Mx[1] = ifelse(x < C[1], 1, 0)
  Mx[J] = ifelse(x >= C[J-1], 1, 0)
  for(j in 2:(J-1))
    Mx[j] = ifelse((C[j-1] <= x) & (x < C[j]), 1, 0)
  return (Mx)
}
## categorized covariate function of x: x is a vector.
m <- function(x){
  n = length(x)
  # Categorize W
  cx = matrix(0, nrow = n,ncol = J)
  for(i in 1:n){cx[i, ] = M(x[i])}
  return(cx)}
## the row sum of misclassification matrix is 1.
fn_norm <- function(A) {
  A[A < 0] = 0
  A = A / rowSums(A)
  return(A)
}
## compute the n-th power of A
fn_power <- function(A,n){
  a = eigen(A)$values
  Sigma = diag(a^n)
  P = eigen(A)$vectors
  out = P%*%Sigma%*%solve(P)
  return(out)
}
## fn_lambda is used in fn_simfex
fn_lambda <- function(data,A,p,lambda,theta){
  p = t(A)%*%p
  A = fn_power(A,lambda)
  A = fn_norm(A)
  theta_lam = c()
  for(i in 1:J){theta_lam[i] =  sum(A[,i]*p*theta)/sum(A[,i]*p)}
  return(theta_lam)
}
## function for our simfex method
fn_simfex <- function(w,y,B.boot){
  theta_naive = theta_simfex = rep(NA,J+1)
  se.naive = b.se.naive = b.se.simfex = rep(NA,J+1)
  b.theta.naive = b.theta.simfex = matrix(NA,B.boot,J+1)
  # 1.estimate nuisance parameters, misclassification matrix A_box and probability vector p------------
  ## define cut points: given
  ## estimate lambda in Box-Cox transformation
  b <- boxcox(lm(w ~ 1))
  lam_hat <- b$x[which.max(b$y)]  
  ## estimate mean_box_x, s2_u, s2_box_x
  box_W_int = matrix(0,n,r)
  for(j in 1:r){
    box_W_int[,j] = Box(W_int[,j],lam_hat)
  }
  box_W_int = na.omit(box_W_int)
  row_mean_box_w = apply(box_W_int, 1, mean) 
  mean_box_x_hat = mean(row_mean_box_w) 
  s2_box_w = apply(box_W_int, 1, var) 
  s2_u_hat = mean(s2_box_w) 
  s2_box_x_hat = max(mean((row_mean_box_w - mean_box_x_hat) ^ 2) - s2_u_hat/r, 
                     0.2*(mean((row_mean_box_w - mean_box_x_hat)^2))) 
  ## estimate misclassification matrix A_box and probability vector p
  a = rep(0,(J+1)); a[1] = min(w); a[2:J] = C; a[J+1] = max(w)
  Ap = fn_box(mean_box_x_hat,s2_box_x_hat,s2_u_hat,lam_hat,a)
  A_box = Ap$A_box 
  p_hat = Ap$p_hat  
  
  # 2.estimate theta_J-theta1 using naive, simfex methods------------
  ## input
  w0 = w # observations with measurement error, a n-dim vector
  W0 = W_int # replicate observations with measurement error, a n*r matrix
  y0 = y # response, a n-dim vector
  A = A_box # a misclassification matrix with J rows and J columns
  p = p_hat # a J-dim probability vector
  B.boot = B.boot 
  ## naive estimator
  ### Run standard logistic regression using glm (no intercept)
  thetaw_out = glm(y ~ m(w) - 1, family = binomial(link = "logit"))
  out = summary(thetaw_out)$coef[1:J]
  out = c(out,out[J]-out[1])
  theta_naive = out
  ## simfex estimator
  lambda = c(0.5,1,1.5,2)
  theta = matrix(NA,length(lambda),J)
  for(l in 1:length(lambda)){
    a_lam = fn_lambda(m(w),A,p,lambda[l],H_fn(theta_naive[-(J+1)]))
    theta[l,] = a_lam
  }
  par(mfrow=c(3,2))
  out = c()
  for(j in 1:J){
    plot(lambda,theta[,j],main=paste("simfex:theta_",j,sep=""))
    ### extrapolation function(quadratic)
    fit = lm(theta[,j] ~ lambda + I(lambda^2))
    a = as.vector(fit$coefficients)
    out[j] = a[1]-a[2]+a[3]
  }
  out = H_fn_inv(out)
  out = c(out,out[J]-out[1])
  theta_simfex = out
  
  ## bootstrap to get variance of naive and simfex----------------
  for(b in 1:B.boot){
    b.data = matrix(0,n,1+r)
    b.data[,1] = y0
    b.data[,2:(1+r)] = W0
    index = sample(1:n,size = n,replace = TRUE)
    b.data = b.data[index,]
    y = b.data[,1]
    W_int = as.matrix(b.data[,2:(1+r)])
    w = W_int[,1]
    ## reestimate nuisance parameters, misclassification matrix A_box and probability vector p
    b.boxcox <- boxcox(lm(w ~ 1))
    lam_hat <- b.boxcox$x[which.max(b.boxcox$y)]
    box_W_int = matrix(0,n,r)
    for(j in 1:r){
      box_W_int[,j] = Box(W_int[,j],lam_hat)
    }
    box_W_int = na.omit(box_W_int)
    row_mean_box_w = apply(box_W_int, 1, mean)
    mean_box_x_hat = mean(row_mean_box_w)
    s2_box_w = apply(box_W_int, 1, var)
    s2_u_hat = mean(s2_box_w)
    s2_box_x_hat = max(mean((row_mean_box_w - mean_box_x_hat) ^ 2) - s2_u_hat/r,
                       0.2*(mean((row_mean_box_w - mean_box_x_hat)^2)))
    a = rep(0,(J+1)); a[1] = min(w); a[2:J] = C; a[J+1] = max(w)
    Ap = fn_box(mean_box_x_hat,s2_box_x_hat,s2_u_hat,lam_hat,a)
    A = Ap$A_box
    p = Ap$p_hat
    ## naive --------------------------------
    thetaw_out = glm(y ~ m(w) - 1, family = binomial(link = "logit"))
    out = summary(thetaw_out)$coef[1:J]
    out = c(out,out[J]-out[1]) 
    b.theta.naive[b,]= out
    ## simfex ------------------------------------------------------
    lambda = c(0.5,1,1.5,2)
    theta = matrix(NA,length(lambda),J)
    for(l in 1:length(lambda)){
      a_lam = fn_lambda(m(w),A,p,lambda[l],H_fn(b.theta.naive[b,][-(J+1)]))
      theta[l,] = a_lam}
    par(mfrow=c(3,2))
    out = c()
    for(j in 1:J){
      plot(lambda,theta[,j],main=paste("simfex:theta_",j,sep=""))
      ### extrapolation function(quadratic)
      fit = lm(theta[,j] ~ lambda + I(lambda^2))
      a = as.vector(fit$coefficients)
      out[j] = a[1]-a[2]+a[3]
    }
    out = H_fn_inv(out)
    out = c(out,out[J]-out[1])
    b.theta.simfex[b,] = out
  }
  b.se.simfex = apply(b.theta.simfex, 2, sd)

  # 3. output-----------------------------------------------------------------
  ## estimate, se, p.value (H_0: theta_J-theta_1 = 0)
  theta_simfex = theta_simfex[J+1] 
  se.simfex = b.se.simfex[J+1]
  p.simfex = pnorm(abs(theta_simfex/se.simfex), lower.tail = F)*2
  result = matrix(round(c(theta_simfex, se.simfex, p.simfex),3),1,3)
  rownames(result) = "theta_J-theta_1"
  colnames(result) = c("estimate", "se", "p.value")
  return(result)
}

