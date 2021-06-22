library(class)

S_func <- function(x, a) {  # Soft Thresholding Operator
  
  
  return( pmax(abs(x) - a,0) * sign(x))
}

#library(RGCCA)# soft threshold
#_func=soft.threshold
Log<-function(y,pr){
  
  if(isTRUE(y==0)==T){a<-0
  
  } else{
    a<- y*log((y/pr))
  }
 
  return(a)
}



errfun.binomial=function(y,yhat,w=rep(1,length(y))){
  prob_min = 1e-05
  prob_max = 1 - prob_min
  predmat = pmin(pmax(yhat, prob_min), prob_max)
  -w*(y*log(predmat)+(1-y)*log(1-predmat))
}






######
library(pracma)

reg<-function(r,Z){
  K=ncol(Z)
  my_one<-matrix(1,nrow(Z))
  my_w=data.frame(Z,my_one)
  my_w<-as.matrix(my_w)
  my_inv<-pinv(t(my_w)%*%my_w)
  my_res<-my_inv%*%(t(my_w)%*%r)
  # new<- lm(r~1,na.action=na.exclude)
  beta0<-matrix(my_res[(K+1)])
  # new1<- lm(r~Z,singular.ok = TRUE)
  
  theta0<- matrix(my_res[c(1:(K))])
  
  return(list(beta0,theta0))
}




error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  #range(upper, lower)
}



twonorm<-function(x){
  sqrt(sum(x*x))
  
}


objective<-function(r,beta,theta,alpha,lambda){
  
  p=length(beta)
  norm_1=   lapply(seq_len(p),
                     function(g)(  lambda*(1-alpha)*(norm(matrix(c(beta[g],theta[g,])),type = "F") +norm(matrix(c(theta[g,])),type = "F") ) + lambda*alpha* sum(abs(theta[g,]))      ))
  

  objective_1 <- sum(r^2) / (2 * N) +
    sum(unlist(norm_1))
  return(objective_1)
}






gradient_j<-function(beta,theta,U,U2,U3,y,X,W,r,alpha,lambda,K,j,N){
  
  theta_j<- matrix(as.numeric(theta[j,]),nrow = 1)
  
  dJ_dtheta_j<-matrix(0,nrow = 1,ncol = K)
  R<-c()
  
  if (isTRUE(any(c((beta[j]), theta_j)!=0))==T) {
    u<-as.numeric(beta[j])/norm(matrix(c(beta[j], theta_j)),type= "F")
    U<-c(U,u)
    
  } else{
    if(length(U)>=1){
      for (z in 1:length(U)) {
        if (norm(matrix(U[z]),type = "F")<=1) {
          R<-c(R,U[z])
          next(z)
        } else{
          next(z)
        }
        
        
      }
    }else{R=c(R,0)}
    u<-sample (c(R), size=1)
    U<-c(U,u)
  }
  
  
  
  
  dJ_dbeta_j=-(t(matrix(X[,j] )) %*% matrix((matrix(r)) ) )/N+(1-alpha)*lambda*u
  #DJ_BETA[j]<-dJ_dbeta_j
  
  
  
  
  
  #for (n in 1:K) {
  
  R<-c()
  if (isTRUE(any(c((beta[j]), theta_j)!=0))==T) {
    u2<-as.numeric(theta_j)/norm(matrix(c(as.numeric(beta[j]), theta_j)),type= "F")
    U2<-c(U2,u2)
    
  } else{
    if(length(U2)>=1){
      for (z in 1:length(U2)) {
        if (norm(matrix(U2[z]),type = "F")<=1) {
          R<-c(R,U2[z])
          next(z)
        } else{
          next(z)
        }
        
        
      }
      
    }else{R=c(R,0)}
    u2<-sample (c(R), size=1)
    U2<-c(U2,u2)
  }
  
  
  d<-c()
  if (isTRUE(any(c(theta_j)!=0))==T) {
    u3<-as.numeric(theta_j)/norm(matrix(c(theta_j)),type= "F")
    U3<-c(U3,u3)
    
  } else{
    if(length(U3)>=1){
      for (z in 1:length(U3)) {
        if (norm(matrix(U3[z]),type = "F")<=1) {
          d<-c(d,U3[z])
          next(z)
          
        } else{
          next(z)
        }
        
        
      }
    }else{d=c(d,0)}
    u3<-sample (c(d), size=1)
    U3<-c(U3,u3)
  }
  
  
  
  
  v<-sign(as.numeric(theta_j))
  
  
  
  
  
  
  
  
  dJ_dtheta_j= -((t(data.frame(W[j]) )%*% matrix(matrix((r) ) ) ))/N+(1-alpha)*lambda*(u2+u3)+alpha*lambda*v
  #}
  
  
  #DJ_THETA[j,]<-dJ_dtheta_j
  
  
  
  
  
  
  L1<- dJ_dbeta_j
  L2<- dJ_dtheta_j
  
  
  return(list(L1,L2,U,U2,U3))
}



compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))
  
  return(xz_term)
  
  
}



concat_beta_theta<-function(beta, theta){
  p=length(beta); K =ncol(theta) 
  my_matrix = matrix(0,p, (K+1))
  my_matrix[, c(1:K)] = theta
  my_matrix[, (K+1)] = beta
  return (my_matrix)
  
}

model<-function(beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  #The pliable lasso model described in the paper
  #y ~ f(X)
  
  #formulated as
  
  #y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )
  
  
  
  intercepts = as.numeric(beta0)+Z%*%(matrix(theta0))
  shared_model = X%*%matrix(beta)
  pliable = compute_pliable(X, Z, theta)
  return( intercepts + shared_model + pliable)
}



model_min_j<-function(beta0, theta0, beta, theta, X, Z, j, W ){
  
  #y ~ f(x) with X_j removed from the model
  
  
  beta[j] = 0.0
  theta[j, ] = 0.0
  return (model(beta0, theta0, beta, theta, X, Z))
  
}



model_j<-function(beta_j, theta_j, x_j, W, j, Z){
  
  
  #r_j ~ beta_j * X_j + W_j @ theta_j
  
  #Only a the residual fit on a single predictor is made
  
  
  # Caching is disabled
  #zz=matrix(compute_w_j(x_j, Z,matrix(theta_j,1)))
  w_j = as.matrix(data.frame(W[j]))
  zz<- w_j%*%theta_j
  return (beta_j * x_j + zz)
}



penalties_min_j<-function(beta_0, theta_0, beta, theta, X, Z, y, W, j,E,n_i,b){
  # Compute the MSE, penalty 1 and penalty 2 when the jth predictor is not in the model.
  N=nrow(X)
  p=ncol(X)
  
  
  #y_hat_min_j = model_min_j(beta_0, theta_0, beta, theta, X, Z, j, W)
  
  n_i_b<-model_min_j(beta_0, theta_0, beta, theta, X, Z, j, W)
  #E1=E
 # E1[,b]<-0
  

  pr_b<-  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
  
  
  
  
  
  
  B<-pr_b*(1-pr_b)
  
  
  B[B==0]=as.numeric(10e-9)
  
  
  y1=1*(y==b)
  
  M<-(y1-pr_b)
  
  A<-n_i_b+M/B
  rv<-B*(A-n_i_b)
  
  
  
  mse = (1/(2*N))/ sum((rv)^2)
  
  coef_matrix = concat_beta_theta(beta, theta)
  
  # Ignore the jth modifier from the model
  coef_matrix[j, ] = 0.0
  theta[j, ] = 0.0
  
  # Compute penalties
  penalty_1=0; penalty_2 = 0.0
  
  for (j in 1: p){
    penalty_1= penalty_1 + twonorm(matrix(coef_matrix[j, ]))
    penalty_2=penalty_2 + twonorm(matrix(theta[j, ]))
  }
  penalty_3 = sum(abs(theta))
  
  
  return (list(mse, penalty_1, penalty_2, penalty_3))
  
}



objective_j<-function(beta0,theta0,beta,theta,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j,b,n_i,E,B){
  beta_j<-beta[j]
  theta_j<-theta[j,]
  theta_transpose<-t(theta_j)
  K=length(theta_j)
  # print(j)
  
  precomputed_penalties_minus_j = penalties_min_j(beta0, theta0, beta, theta, X, Z, y, W, j,E,n_i,b)
  # 
  # 
  mse<-as.numeric(unlist(precomputed_penalties_minus_j[1]))
  penalty_1<-as.numeric(unlist(precomputed_penalties_minus_j[2]))
  penalty_2<-as.numeric(unlist(precomputed_penalties_minus_j[3]))
  penalty_3<-as.numeric(unlist(precomputed_penalties_minus_j[4]))
  
  
  #r_min_j = ((model_min_j(beta0, theta0, beta, theta, X, Z, j, W )))
  
  # print(theta_j)
  # zz<- as.matrix(data.frame(W[j]))%*%theta_j
  #print(zz)
  r_hat = model_j(beta[j], theta[j,], X[, j], W, j, Z)
  mse_1 = (1 / (2*N)) * ( sum(  (r_min_j - B*r_hat )^2 )) +mse
  
  
  # Penalty 1
  coef_vector = matrix(0,(K+1))
  coef_vector[c(1:K)] = as.numeric(theta_j)
  coef_vector[(K+1)] = beta_j
  penalty_1 =penalty_1+ twonorm(matrix(coef_vector) )
  
  # Penalty 2
  penalty_2 =penalty_2+ twonorm(matrix(theta_j))
  
  # Penalty 3
  penalty_3 =penalty_3+ sum(abs(theta_j))
  
  objective_l = mse_1 + (1-alpha) * lambda * (penalty_1 + penalty_2) + alpha * lambda * penalty_3
  
  
  
  
  
  
  
  
  # norm_1_l=     lambda*(1-alpha)*(norm(matrix(c(beta_j,theta_j)),type = "F") +norm(matrix(c(theta_j)),type = "F") ) + lambda*alpha* sum(abs(theta_j))      
  
  #objective_l=apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)
  
  
  
  return( objective_l)
}



# 
# res_min_j<-function(beta_0, theta_0, beta, theta, X, Z, y, W,b,n_i,E,pr,j){
#   # Compute the MSE, penalty 1 and penalty 2 when the jth predictor is not in the model.
#   N=nrow(X)
#   p=ncol(X)
#   
#   
#   y_hat_min_j = model_min_j(beta_0, theta_0, beta, theta, X, Z, j, W)
#   n_i_b<-model_min_j(beta_0, theta_0, beta, theta, X, Z, j, W)
#   E1=E
#   E1[,b]<-0
#   
#   pr[b]<-list(matrix(as.numeric(1/( exp(-n_i_b)*rowSums(E1)+1),N)))
#   pr_b<-matrix(unlist(pr[b]),N)
#   
#   
#   
#   
#   
#   
#   B<-pr_b*(1-pr_b)
#   
#   
#   B[B==0]=as.numeric(10e-9)
#   
#   
#   y1=1*(y==b)
#   
#   M<-(y1-pr_b)
#   
#   A<-n_i_b+M/B
#   rv<-B*(A-n_i_b)
#   
#   
#   
#   return (rv)
#   
# }
# 


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
####################TRAIN TRUE WORK#################################

for_v=10;sv=0;fq=50;st=50;mv=20;ms=50;tol=1e-3



plasso_fit1 <- function(y, X, Z, nlambda, alpha,new_t,my_mbeta,number, intercept,step,maxgrid,tol,run,lambda_min,my_lambda=NULL,tt=NULL,for_v,sv,fq,st,mv,ms,cv_run,max_iter=100){
  
 # xbar=colMeans(X)
 #zbar=colMeans(Z)
  
 #X=scale(X,xbar,F)
 # Z=scale(Z,zbar,F)
 # print(Z)
  orig.X=X
  #y=as.vector(y)
  N <- length(y)
  
  p <- ncol(X)
  K <- ncol(Z)
  #X=matrix(as.numeric(X),N,p)
 # Z=matrix(as.numeric(Z),N,K)
  #if (p<N) {rat<-1e-5
  
  # } else {rat<-1e-2} 
  
  rat=lambda_min
  
  
  tolerance=tol
  
  
  W <- lapply(seq_len(p),
              function(j) (matrix(X[, j], nrow = N, ncol = K) * Z))
  
  
  
  
  
  #for (n in 1:max(y)) {
  # BETA0[n]<-data.frame(matrix(0)); BETA[n]<-data.frame(matrix(0)); THETA0[n]<-data.frame(matrix(0)); THETA[n]<-data.frame(matrix(0)); Y_hat[n]<-data.frame(matrix(0)); Lambda[n]<-data.frame(lambda); non_zero[n]<-data.frame(matrix(0, nrow= nlambda))
  # }
  
  
  # BETA0<-matrix(0,nrow = (nlambda)); BETA<-matrix(0,nrow = (nlambda)); THETA0<-matrix(0,nrow = (nlambda)); THETA<-matrix(0,nrow =( nlambda)); 
   
   
   Y_hat<-matrix(0, nrow = (nlambda)); Lambda<-matrix(0, nrow = (nlambda)); non_zero<-matrix(0, nrow= (nlambda),ncol = max(y));non_zero_theta<-matrix(0, nrow= (nlambda),ncol = max(y));DEV<-matrix(0, nrow = (nlambda),ncol = 1);DEV1<-matrix(0, nrow = (nlambda),ncol = 1)
  # 
  
  
  BETA0<-lapply(seq_len(max(y)),
                 function(j)(matrix(0,nrow = (nlambda))))
  
  BETA<-lapply(seq_len(max(y)),
                function(j)(matrix(0,nrow=p,ncol=nlambda)))
  
  THETA0<-lapply(seq_len(max(y)),
                  function(j)(matrix(0,nrow=K,ncol=nlambda)))
  
  THETA<-lapply(seq_len(max(y)),
                 function(j)(array(0,c(p,K,nlambda))))
  
  #y_Hat <-  X%*%matrix(beta) + XZ_termnggb
  
  
  
  
  #beta01<-beta; theta01<-theta
  
  
  
  
  
  # BETA0[1]<-list(beta0); BETA[1]<-list(matrix(0,p,max(y))); THETA0[1]<-list(theta0); THETA[1]<-list(theta)
  
  
  
  # while (lambda>=Lambda_min ){
  # 
  BETA01<-lapply(seq_len(max(y)),
                 function(j)(0))
  old_BETA01<-lapply(seq_len(max(y)),
                     function(j)(0))
  BETA1<-lapply(seq_len(max(y)),
                function(j)(matrix(0,nrow = 1,ncol = p)))
  old_BETA1<-lapply(seq_len(max(y)),
                    function(j)(matrix(0,nrow = 1,ncol = p)))
  #
  THETA01<-lapply(seq_len(max(y)),
                  function(j)(matrix(0,nrow = 1,ncol = K)))
  old_THETA01<-lapply(seq_len(max(y)),
                      function(j)(matrix(0,nrow = 1,ncol = K)))
  THETA1<-lapply(seq_len(max(y)),
                 function(j)(matrix(0,nrow = p,ncol = K)))
  old_THETA1<-lapply(seq_len(max(y)),
                     function(j)(matrix(0,nrow = p,ncol = K)))
  Y_hat1<-lapply(seq_len(max(y)),
                 function(j)(matrix(0,nrow = N)))
  
  XZ_term<-lapply(seq_len(max(y)),
                  function(j)(matrix(0,nrow = N)))
  n_i<-lapply(seq_len(max(y)),
              function(j)(matrix(0,nrow = N)))
  pr<-lapply(seq_len(max(y)),
             function(j)(matrix(0,nrow = N)))
  n_il<-lapply(seq_len(max(y)),
               function(j)(matrix(0,nrow = N)))
  n_ir<-lapply(seq_len(max(y)),
               function(j)(matrix(0,nrow = N)))
  prl<-lapply(seq_len(max(y)),
              function(j)(matrix(0,nrow = N)))
  prr<-lapply(seq_len(max(y)),
              function(j)(matrix(0,nrow = N)))
  
  
  for (x in 1:max(y)) {
    #beta0<-unlist(BETA01[x])
    # beta <- matrix(unlist(BETA1[x]),p)
    # my_X<-scale(X,-xbar,FALSE) ; my_Z<-scale(Z,-zbar,FALSE)
    my_X<-X ; my_Z<-Z
    #my_X=matrix(as.numeric(my_X),N,p)
    #my_Z=matrix(as.numeric(my_Z),N,K)
    
    
    beta0<-BETA0[[x]][1]
    beta <- BETA[[x]][,1]
    # print(beta)
    
    theta <-as.matrix(THETA[[x]][,,1]) ; theta_transpose <- t(theta)
    
    
    theta0 <-THETA0[[x]][,1]
    
    
    
    
    
    
    
    
    # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
    
    n_i[x]<-list(model(beta0, theta0, beta, theta, X, Z))
    
  }
  E<-matrix(0,N,max(y))
  for (x in 1:max(y)) {
    E[,x]<-exp(unlist(n_i[x]))
  }
  my_pr<-matrix(0,N,(max(y)-1))
  for (x in 1:(max(y))) {
    n_i_b<-matrix(unlist(n_i[x]),N)
    #E_1<-E
    #E_1[,x]<-0
    pr[x]<-list(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
    #my_pr[,x]<-(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
  }
 # pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
  
  
  
  
  
  
  
  if(is.null(my_lambda)){
    
    O_1<-matrix(0,max(y),1)
    
    #active_set1<-matrix(0,p,max(y))
    
   # strong_set<-matrix(0,p,max(y))
    #strong_set1<-matrix(0,p,max(y))
    
    
    for(b in 1:max(y)) {
      n_i_b<-matrix(unlist(n_i[b]),N)
      #beta0<-matrix(unlist(BETA01[b]))
      pr_b<-matrix(unlist(pr[b]),N)
      
      
      
      B<-pr_b*(1-pr_b)
      B[B==0]=as.numeric(10e-9)
      # O<-B
      
      
      y1=1*(y==b)
      M<-(y1-pr_b)
      
     A<-n_i_b+M/B
     # rv<-B*(A-n_i_b)
     r<-B*(A-n_i_b)
     
      #
     # r=y1
      ## Define_Parameters
      #beta01<-beta0
      O_1[b]=max(abs(t(X)%*%r)/length(r))/(1-alpha)
      #print(O)
     
      
      
      
     
      
    }
    print(O_1)
   # lambda_to_use<-matrix(0,nlambda,max(y))
    #for (i in 1:max(y)) {
    
    
    # lambda_max_to_select<-matrix(0,max(y))
    
    
    #for (j in 1:max(y)) {
    lambda_max_to_select<-   max(O_1)
    #}
    lambda<-max(lambda_max_to_select)
    print(lambda)
    big_lambda<-lambda
    #lambda_i<-matrix(0,maxgrid,max(y))
    
    Lambda_min<- rat*big_lambda
    print(Lambda_min)
   # lambda_i<-matrix(0,maxgrid)
    #for (g in 2:(maxgrid-1)) {
    #lambda_i[g]<- big_lambda*(Lambda_min/big_lambda)^(g/maxgrid)
    #}
    lambda_i<- exp(seq(log(big_lambda),log(big_lambda*rat),length=maxgrid))
    lambda_i[1]<-big_lambda;lambda_i[maxgrid]<-Lambda_min
        
     # lambda_i<-exp(seq(log(big_lambda),log(big_lambda*rat), (log(big_lambda*rat) - log(big_lambda))/(maxgrid-1)))
print(lambda_i)
    #lambda_i[1]<-big_lambda ;lambda_i[maxgrid]<-Lambda_min
    #lambda_i[1]<-lambda_max_to_select;lambda_i[maxgrid]<-Lambda_min
    # lambda_to_use[,i]<-lambda_i
    # }
  }else{
    lambda_i<-my_lambda
  }
  
  big_t<-matrix(0,nlambda,max(y))
  Mbeta<-matrix(0,nlambda)
  #lambda<-lambda_i[1]
  #t=new_t
  mbeta=my_mbeta
  
  #v_beta<-0
  #v_theta<-matrix(0,1,K)
  #old_beta_hat<-matrix(0,p,1)
  #old_theta_hat<-matrix(0,p,K)
  #nes_lambda<-matrix(0,nlambda+1)
  
  #for (i in 1:nlambda) {
  #  nes_lambda[1]<-0
  #  nes_lambda[i+1]<-(1+sqrt(1+4*(nes_lambda[i])^2))/2
  #}
  dev_percent<-matrix(0,nlambda)
  sec_active1<-matrix(0,p,max(y))
  sec_active2<-matrix(0,p,max(y))
  
  ACTIVE<-matrix(0,nlambda,max(y))
  ACTIVE1<-matrix(0,nlambda,max(y))
  my_v<-sv
  my_V<-for_v
  
  active_set1<-matrix(0,p,max(y))
  active_set2<-matrix(0,p,max(y))
 
  strong_set<-matrix(0,p,max(y))
  strong_set1<-matrix(0,p,max(y))
  
  my_q=1
  my_ok<-0
  q=1
  
  U<-matrix(0); U2<-matrix(0); U3<-matrix(0)
  while (isTRUE(my_ok==0)==T | isTRUE(q<=nlambda)==T ) {
  #for (i in 1:nlambda ) {
   
    my_v<-my_v
    my_V<-my_V
   # print(c(my_v,my_V))
    
    if(q<=1){
      
      i=q
      
      lambda<-lambda_i[i]
      
      
      
     
      
      
      
      for (b in 1:max(y)) {
        # lambda_i<- lambda_to_use[,b]
      #  for (it in 1:max_it) {
          
       
        
        
        
        beta0=BETA0[[b]][q]
        beta = BETA[[b]][,q]
        # print(beta)
        v_beta=beta
        theta =(THETA[[b]][,,q]) ; theta_transpose <- t(theta)
        v_theta=theta;theta_transpose1 = t(v_theta)
        
        theta0 = THETA0[[b]][,q]
        theta01=theta0
        
        beta01=beta0
        beta1 = beta
        theta1 = theta ; theta_transpose1 = t(theta1)
        theta01 = theta
        norm1<-matrix(0,p,1);norm3<-matrix(0,p,1)
        
        
        for( iii in 1: max_iter){
          # iter_prev_score = objective(
          #   beta0, theta0, beta, theta,
          #   X, Z, y,
          #   alpha, lam, W
          # )
          
          
          
          # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
          n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
         
        
        n_i_b<-matrix(unlist(n_i[b]),N)
       # E1=E
        #E1[,b]<-0
        E[,b]<-exp(unlist(n_i[b]))
          pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
        pr_b<-matrix(unlist(pr[b]),N)
        
        
        
        
        
        
        B<-pr_b*(1-pr_b)
        
        
         B[B==0]=as.numeric(10e-9)
      
        
        y1=1*(y==b)
        
          M<-(y1-pr_b)
        
        A<-n_i_b+M/B
        rv<-B*(A-n_i_b)
        #objective<-function(r,beta,theta,alpha,lambda)
        iter_prev_score= objective(rv,beta,theta,alpha,lambda)
          #iter_prev_score = objective(rv,beta,theta,alpha,lambda )
        
        
        #r<-(A-n_i_b)
        v=reg(rv,Z)
        
        if (b==max(y)) {
          my_beta0<-matrix(0,(max(y)-1))
          for (l in 1:(max(y)-1)) {
            my_beta0[l]<-BETA0[[l]][i]
          }
          beta0<-1-sum(my_beta0)  
          theta0<-matrix(unlist(v[2]))
        } else{
          
          beta0<-matrix(unlist(v[1]))
          theta0<-matrix(unlist(v[2]))
          }
        
        
        
      # theta03<- matrix(0,N,K)
      # xz_theta <- lapply(seq_len(p),
       #                  function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% theta_transpose[, j])
        
     # XZ_term <- Reduce(f = '+', x = xz_theta)
        
        
     #  for (k in 1:K) {
      #  theta03[,k]<- (Z[,k]*(B)*Z[,k])
      # }
     #  theta03<- colSums(theta03)
     #   theta03<- (matrix(theta03))^-1
      #  theta02<-t(Z)%*%((B)*(A-(as.numeric(beta0)+X%*%matrix(beta) + XZ_term)))
       # theta0<- theta03*theta02
        
        
        
        
        #ACTIVE[i,b]<-sum(strong_set)
        
        
        #######################line search###########   
        for (j in 1:p) {
          
          if(is.null(tt)){
            t=NULL
            
          } else{t=tt[q]}
          
         
          # Strong screening rule (https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
          if(i<=1){checkk=2*lambda-lambda_i[1] }else{checkk=( 2*lambda-lambda_i[i-1] )}
          
          
          
          # print(c(i,as.numeric(abs( (t(matrix(X[,j]*c(y) ))%*% matrix((r1) ) )  )/(N*alpha) ), checkk ))
          
          if(as.numeric(abs( (t(matrix( (X[,j]) ) )%*% matrix(((rv))   ) )  )/(N*alpha) )< checkk  ){
            next(j)
          } else { 
            #w_j = ( data.frame(W[j]) )
           res_j<- rv+B*model_j(beta[j], theta[j,], X[,j], W, j, Z) 
           
           cond1<- as.numeric(abs((t(matrix(X[,j]))%*% matrix( res_j ))/N))
            
            
            cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[j]))%*%( res_j ) )/N, alpha*(lambda))), type="F"))
            
            
            
           # cond1<- as.numeric(abs((t(matrix(X[,j]))%*% matrix(y1-pr_b))/N))
            
            
            #  cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[j]))%*%(y1-pr_b))/N, alpha*(lambda))), type="F")) 
            
            
            # cond1<- as.numeric(abs((t(matrix(X[,v]))%*% matrix(r_j))/N))
            
            
            #  cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[v]))%*%(r_j))/N, alpha*(lambda))), type="F")) 
            
            #lambda_i[0]<-0
            #if(isTRUE(active_set1[v,b]==1)==T){
            # strong_set[v]<-1  
            # next(v)
            # } else
            if (cond1<= (1-alpha)*(lambda)){
              
              strong_set[j,b]<-0
            } else{
              strong_set[j,b]<-1
            }
            # print(c(cond1,(1-alpha)*(lambda)))
            
            if (cond2<= 2*(1-alpha)*(lambda)){
              
              strong_set1[j,b]<- 0
            } else{
              strong_set1[j,b]<-1
            }
            
               
             # print(c(cond1,(1-alpha)*(lambda),cond2,2*(1-alpha)*(lambda)))
            if(strong_set[j,b]==0 & strong_set1[j,b]==0){
              
              #beta[j]<-0;theta[j,]<-0
              next(j)
              
            }else  {  
            
            
            
            
           # beta1_j<- (N/(t(B)%*%(matrix(X[,j]))^2))*  S_func((t(B*matrix(X[,j]))%*%matrix(r_j))/N,(1-alpha)*lambda)
            
            
            
            
           # cond3<- norm(matrix(S_func((t(data.frame(W[j]))%*%((r_j)-(B)*matrix(X[,j])*as.numeric(beta1_j))/N),alpha*(lambda))),type="F")
           
            
          
            beta1_j<- (N/sum(B*(matrix(X[,j])*matrix(X[,j]))))* S_func((t(matrix(X[,j]))%*%matrix(res_j ))/N,(1-alpha)*lambda)
           # #  
           # # 
           #print(beta1_j)
           # #  
           # #       
           #  }
          
            
           
            
            
            
            
              
              
              cond3<- norm(matrix(S_func(t(data.frame(W[j]))%*%( ( res_j )-(B)*matrix(X[,j])*as.numeric(beta1_j))/N,alpha*lambda)),type="F")
              
            if (isTRUE(cond3<= 2*(1-alpha)*lambda)==T) {
              beta[j]<-as.numeric(beta1_j)
             # theta[j,]<-as.numeric(0)
              active_set1[j,b]<-1
              # v_beta[j]<-unlist(value1[[1]]); v_theta[j,]<- unlist(value1[[2]]); beta[j]<-unlist(value1[[3]]);theta[j,]<-unlist(value1[[4]])
              
              
              next(j)
              
              
            }  else {
              
              
          
              
             
              
              # L=mbeta; mu=intercept
              #h=t; my_b= (1-nes_lambda[i])/nes_lambda[i+1]
             # quadratic=function(L1,L2, beta_j ,theta_j, alpha,lambda,t,mbeta,v_beta,v_theta,my_n,beta0,theta0,j,W,X,Z,y1,y,N)
             
              
              value1<-   quadratic(beta,theta, alpha,lambda,beta0,theta0,j,b,W,X,Z,y,N,n_i,xbar,zbar,t=t,r_min_j = res_j)
              #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
              
              beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]])
              
              # v_beta[j]<-unlist(value1[[1]]); v_theta[j,]<- unlist(value1[[2]]); beta[j]<-unlist(value1[[3]]);theta[j,]<-unlist(value1[[4]])
              active_set1[j,b]<-1
              active_set2[j,b]<-1
              
              next(j)} #theta
            
          }
          
          }
          
        }#J
        
      
        
        
        # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
        n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
        
        
        n_i_b<-matrix(unlist(n_i[b]),N)
        E[,b]<-exp(unlist(n_i[b]))
        pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
        pr_b<-matrix(unlist(pr[b]),N)
        
        
        
        
        
        
        B<-pr_b*(1-pr_b)
        
        
        B[B==0]=as.numeric(10e-9)
        
        
        y1=1*(y==b)
        
        M<-(y1-pr_b)
        
        A<-n_i_b+M/B
        rv<-B*(A-n_i_b)
        tolerance=tol
        if(i>new_t){
          # iter_current_score = objective(
          #   beta0, theta0, beta, theta,
          #   X, Z, y=orig_y,
          #   alpha, lam, W
          # )
          
          iter_current_score= objective(rv,beta,theta,alpha,lambda)
          #print(c(iter_prev_score,iter_current_score,beta))
          if(isTRUE(abs(iter_prev_score - iter_current_score) <tolerance)==T ){
            # print(c(iter_prev_score,iter_current_score))
            break  # Converged on lam_i
            
          }else{next(iii)}
          
        }else{break}
        #objective<-function(r,beta,theta,alpha,lambda)
       
        
        
        
      }#max it
        
        BETA01[b]<-list(beta0) 
        
        BETA1[b]<-list(beta)
        
        THETA01[b]<-list(matrix(theta0,1,K))
        THETA1[b]<-list(as.matrix(theta,p,K))
        
        BETA0[[b]][q]<-beta0
        THETA0[[b]][,q]<-theta0
        BETA[[b]][,q]<-beta
       # print(beta1)
       # print(BETA[[b]][,i])
        THETA[[b]][,,q]<-theta
        
        
        
        #main_norm[b]<- (1-alpha)*lambda*sum(norm1)+ alpha*lambda*sum(norm3)
        
        next(b)}
      
      
      
      # 
       m_for_beta<-matrix(0,p,max(y))
       for (l in 1:max(y)) {
         m_for_beta[,l]<-BETA[[l]][,q]
       }
      # 
       for (l in 1:nrow(m_for_beta)) {
         med<-mean(m_for_beta[l,])
         m_for_beta[l,]<-m_for_beta[l,]-med
       }
      # 
       #BETA<-lapply(seq_len(max(y)),
         #            function(j)(as.matrix(m_for_beta[,j],nrow = 1,ncol = p)))
      
      
       for (l in 1:max(y)) {
         BETA[[l]][,q]<-m_for_beta[,l]
       }
      
      
      for (x in 1:max(y)) {
        #beta0<-unlist(BETA01[x])
        # beta <- matrix(unlist(BETA1[x]),p)
       # my_X<-scale(X,-xbar,FALSE) ; my_Z<-scale(Z,-zbar,FALSE)
       my_X<-X ; my_Z<-Z
        #my_X=matrix(as.numeric(my_X),N,p)
        #my_Z=matrix(as.numeric(my_Z),N,K)
        
       
       beta0<-BETA0[[x]][q]
       beta <- BETA[[x]][,q]
       # print(beta)
      
       theta <-as.matrix(THETA[[x]][,,q]) ; theta_transpose <- t(theta)
      
       
       theta0 <-THETA0[[x]][,q]
       
       
       
      
          
          
        
        
        # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
        n_i[x]<-list(model(beta0, theta0, beta, theta, X, Z))
        
      }
      E<-matrix(1,N,max(y))
      for (x in 1:max(y)) {
        E[,x]<-exp(unlist(n_i[x]))
        
        
      }
      my_pr<-matrix(0,N,(max(y)-1))
      for (x in 1:(max(y))) {
        n_i_b<-matrix(unlist(n_i[x]),N)
       # E_1<-E
       # E_1[,x]<-0
        pr[x]<-list(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
       # my_pr[,x]<-(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
      }
      
      
      #pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
      
      
      #v1_d<-1*(v1==d)
      Dev<-matrix(0,nrow = length(y))
      Dev1<-matrix(0,nrow = length(y))
      for (l in 1:length(Dev)) {
        deviance_y1<-matrix(0,1,max(y))
        deviance_y2<-matrix(0,1,max(y))
        for (d in 1:max(y)) {
          y_d<-1*(y[l]==d)
          # n_i_d<-matrix(unlist(n_i[d]),N)
          prob_d <-matrix(unlist(pr[d]),N)
          # deviance_y1[d]<-y_d*n_i_d[l]
          # deviance_y2[d]<-exp(n_i_d[l])
          deviance_y1[d]<-  Log(y=y_d,pr=prob_d[l]) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
          deviance_y2[d]<-Log(y=y_d,pr=mean(1*(y==d)))
        }
        Dev[l]<-  sum(deviance_y1)
        Dev1[l]<-  sum(deviance_y2)
      }
     # Dev<-matrix(0,nrow = length(y))
     #  Dev1<-matrix(0,nrow = length(y))
     #  for (l in 1:length(Dev)) {
     #   deviance_y1<-matrix(0,1,max(y))
     #  
     # deviance_y2<-matrix(0,1,max(y))
     #  for (d in 1:max(y)) {
     #   y_d<-1*(y==d)
     #  # n_i_d<-matrix(unlist(n_i[d]),N)
     #    prob_d <-matrix(unlist(pr[d]),N)
     #  # deviance_y1[d]<-y_d*n_i_d[l]
     #  # deviance_y2[d]<-exp(n_i_d[l])
     #   deviance_y1[d]<-sum(errfun.binomial(y=y_d[l],yhat=prob_d[l],w=rep(1,length(y_d[l])))) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
     #  
     #    deviance_y2[d]=sum(errfun.binomial(y=y_d[l],mean(y_d),w=rep(1,length(y_d[l]))))
     # }
     # Dev[l]<-  sum(deviance_y1)
     # Dev1[l]<-  sum(deviance_y2)
     #  
     #  
     #  #(y,yhat,w=rep(1,length(y)))
     # }
      DEV1[i]<- ((2)*sum(Dev1))/length(y)
      DEV[i]<- ((2)*sum(Dev))/length(y) #sum(deviance_y1)# #+sum(main_norm)
     #DEV1[i]<-sum(Dev1)
     #  DEV[i]<-sum(Dev)
      # dev_percent[i]<-(sum(deviance_y2)-sum(deviance_y1))/sum(deviance_y2)
      if(is.null(t)){
        t=0
      }else{t=t}
      big_t[i]<- t
      Mbeta[i]<-mbeta
      
      Y_hat[i]<- list(matrix(c(y,unlist(pr)),N,(max(y)+1)))
      
      #y_hat<-list(pr)
      
      
      #Y_hat[q]<-y_hat
      
      for (v in 1:max(y)) {
        
        
        b=matrix(0,p)
        
        
        nzero_beta=  BETA[[v]][,i]
        nzero_theta=unlist(THETA[[v]][,,i])
        # for (m in 1:p) {
        #  if (isTRUE(abs(c[m])>0)==T) {b[m]=1
        
        #  }
        
        #}
        non_zero[i,v]<-(sum(nzero_beta!=0 ) )
        my_ntheta<-matrix(0,p)
        # for (w in 1:p ) {
        #   if(isTRUE(sum( nzero_theta[w,])!=0)==T){
        #     my_ntheta[w]<-1 
        #     
        #   }
        # }
        non_zero_theta[i,v]<-(sum( nzero_theta!=0))
        
      }
      
      #non_zero<- data.frame(non_zero,length(b ) )    
     # BETA[i]<-list(matrix(unlist(BETA1),p,max(y)))
     # BETA0[i]<-list(matrix(unlist(BETA01),1,max(y)))
      #theta01<-theta
      #THETA0[i]<-list(matrix(unlist(THETA01),ncol = max(y)))
      
      #sol=list(beta0,beta,theta0,theta)
     # THETA[i]<-list(THETA1)
      Lambda[i]<-lambda
      #Lambda[i]<-max(lambda_i[i,])
     # print(c(i,t, DEV[i]))
      print(c(q,max(non_zero[i,]) ,big_t[i], DEV[i]))
      
      
      
      my_q<-my_q+1
      q=q+1
      # next(i)
    }else if(isTRUE(q>1)==T & isTRUE(q<=fq)==T){
      
      
      
      i=q
      
      
      lambda<-lambda_i[i]
      
      
      for (b in 1:max(y)) {
        lambda<-lambda_i[q]
        
       # for (it in 1:max_it) {
        
        
        
        beta0<-BETA0[[b]][q-1]
        beta <- BETA[[b]][,q-1]
       # print(beta)
       # print(beta)
        v_beta=BETA[[b]][,q-2]
        theta <-as.matrix(THETA[[b]][,,q-1]) ; theta_transpose <- t(theta)
        v_theta=as.matrix(THETA[[b]][,,q-2]);theta_transpose1 = t(v_theta)
        
        theta0 <- THETA0[[b]][,q-1]
        theta01=theta0
        
        beta01=beta0
        beta1 = beta
        theta1 = theta ; theta_transpose1 = t(theta1)
        theta01 = theta
        norm1<-matrix(0,p,1);norm3<-matrix(0,p,1)
        
        
        for( iii in 1: max_iter){
          # iter_prev_score = objective(
          #   beta0, theta0, beta, theta,
          #   X, Z, y,
          #   alpha, lam, W
          # )
          
          
          
          # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
          n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
          
          
          n_i_b<-matrix(unlist(n_i[b]),N)
          E[,b]<-exp(unlist(n_i[b]))
          pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
          pr_b<-matrix(unlist(pr[b]),N)
          
          
          
          
          
          
          B<-pr_b*(1-pr_b)
          
          
          B[B==0]=as.numeric(10e-9)
          
          
          y1=1*(y==b)
          
          M<-(y1-pr_b)
          
          A<-n_i_b+M/B
          rv<-B*(A-n_i_b)
          #objective<-function(r,beta,theta,alpha,lambda)
          iter_prev_score= objective(rv,beta,theta,alpha,lambda)
          #iter_prev_score = objective(rv,beta,theta,alpha,lambda )
          
       # r<-(A-n_i_b)
        v=reg(rv,Z)
        
        if (b==max(y)) {
          my_beta0<-matrix(0,(max(y)-1))
          for (l in 1:(max(y)-1)) {
            my_beta0[l]<-BETA0[[l]][q-1]
          }
          beta0<-1-sum(my_beta0)  
          theta0<-matrix(unlist(v[2]))
        } else{
          
          beta0<-matrix(unlist(v[1]))
          theta0<-matrix(unlist(v[2]))
        }
        
        
        #r=y1
        
        
        
      
        
        main_norm<- matrix(0,max(y))
        
        
        
        
        
        
      
          
          
         
          #r=y1
          
          
          DJ_BETA<-matrix(0,1)
          DJ_THETA<-matrix(0,1,K)
          #############################################################
          #############################################################
          #############################################################
          #############################################################
          #############################################################    ########## main loop###########
          
          
          
         
          
          for (j in 1:p) {
            if(is.null(tt)){
              t=NULL
              
            } else{t=tt[q]}
            
            # Strong screening rule (https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
            if(i<=1){checkk=2*lambda-lambda_i[1] }else{checkk=( 2*lambda-lambda_i[i-1] )}
            
            
            
            # print(c(i,as.numeric(abs( (t(matrix(X[,j]*c(y) ))%*% matrix((r1) ) )  )/(N*alpha) ), checkk ))
            
            if(as.numeric(abs( (t(matrix( (X[,j]) ) )%*% matrix(((rv))   ) )  )/(N*alpha) )< checkk  ){
              next(j)
            }else  {  
              res_j<- rv+B*model_j(beta[j], theta[j,], X[,j], W, j, Z)
              
              
            cond1<- as.numeric(abs((t(matrix(X[,j]))%*% matrix( res_j ) )/N))
               
              
              cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[j]))%*%( res_j ) )/N, alpha*(lambda))), type="F"))
              
              
              # cond1<- as.numeric(abs((t(matrix(X[,v]))%*% matrix(r_j))/N))
              
              
              #  cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[v]))%*%(r_j))/N, alpha*(lambda))), type="F")) 
              
              #lambda_i[0]<-0
              #if(isTRUE(active_set1[v,b]==1)==T){
              # strong_set[v]<-1  
              # next(v)
              # } else
              if (isTRUE(cond1<= (1-alpha)*(lambda))==T ){
                
                strong_set[j,b]<-0
              } else{
                strong_set[j,b]<-1
              }
              # print(c(cond1,(1-alpha)*(lambda)))
              
              if (isTRUE(cond2<= 2*(1-alpha)*(lambda))==T){
                
                strong_set1[j,b]<- 0
              } else{
                strong_set1[j,b]<-1
              }
                
             # print(c(cond1,(1-alpha)*(lambda),cond2,2*(1-alpha)*(lambda)))
              if(isTRUE(strong_set[j,b]==0)==T & isTRUE(strong_set1[j,b]==0)==T){
                
                #beta[j]<-0;theta[j,]<-0
                next(j)
                
              }else  {  
                
                
                
                
                # beta1_j<- (N/(t(B)%*%(matrix(X[,j]))^2))*  S_func((t(B*matrix(X[,j]))%*%matrix(r_j))/N,(1-alpha)*lambda)
                
                
                
                
                # cond3<- norm(matrix(S_func((t(data.frame(W[j]))%*%((r_j)-(B)*matrix(X[,j])*as.numeric(beta1_j))/N),alpha*(lambda))),type="F")
                
              
               
                beta1_j<- (N/sum(B*(matrix(X[,j])*matrix(X[,j]))))*  S_func((t(matrix(X[,j]))%*%matrix(res_j ))/N,(1-alpha)*lambda)
                #    #  
                #    # 
              # print(beta1_j)
                #    #  
             
                
                
                #print(c(cond3,2*alpha*checkk,beta1_j))
                
             
                
                cond3<- norm(matrix(S_func(t(data.frame(W[j]))%*%( ( res_j )-(B)*matrix(X[,j])*as.numeric(beta1_j))/N,alpha*lambda)),type="F")
                
                #print(c(cond3,2*(1-alpha)*(lambda)))
                if (isTRUE(cond3<= 2*(1-alpha)*(lambda))==T) {
                  beta[j]<-as.numeric(beta1_j)
                  #theta[j,]<-as.numeric(0)
                  active_set1[j,b]<-1
                  #print(active_set1[j,b])
                  # v_beta[j]<-unlist(value1[[1]]); v_theta[j,]<- unlist(value1[[2]]); beta[j]<-unlist(value1[[3]]);theta[j,]<-unlist(value1[[4]])
                  
                  
                  next(j)
                  
                  
                }  else {
                  
                  
         
                  value1<-   quadratic(beta,theta, alpha,lambda,beta0,theta0,j,b,W,X,Z,y,N,n_i,xbar,zbar,t=t,r_min_j = res_j)
                  #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
                  
                  beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]])
                 # print(beta1[j])
                  # v_beta[j]<-unlist(value1[[1]]); v_theta[j,]<- unlist(value1[[2]]); beta[j]<-unlist(value1[[3]]);theta[j,]<-unlist(value1[[4]])
                  active_set1[j,b]<-1
                  active_set2[j,b]<-1
                  
                  next(j)}
              
              }
              
            }
            
          }#j
          
          
        
          
          
          
          # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
          n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
          
          
          n_i_b<-matrix(unlist(n_i[b]),N)
          E[,b]<-exp(unlist(n_i[b]))
          pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
          pr_b<-matrix(unlist(pr[b]),N)
          
          
          
          
          
          
          B<-pr_b*(1-pr_b)
          
          
          B[B==0]=as.numeric(10e-9)
          
          
          y1=1*(y==b)
          
          M<-(y1-pr_b)
          
          A<-n_i_b+M/B
          rv<-B*(A-n_i_b)
          tolerance=tol
          if(i>new_t){
            # iter_current_score = objective(
            #   beta0, theta0, beta, theta,
            #   X, Z, y=orig_y,
            #   alpha, lam, W
            # )
            
            iter_current_score= objective(rv,beta,theta,alpha,lambda)
            #print(c(iter_prev_score,iter_current_score,beta))
            if(isTRUE(abs(iter_prev_score - iter_current_score) <tolerance)==T ){
              # print(c(iter_prev_score,iter_current_score))
              break  # Converged on lam_i
              
            }else{next(iii)}
            
          }else{break}
          #objective<-function(r,beta,theta,alpha,lambda)
          
          
          
          
           }#max it
      
          BETA01[b]<-list(beta0) 
          
          BETA1[b]<-list(beta)
          
          THETA01[b]<-list(matrix(theta0,1,K))
          THETA1[b]<-list(as.matrix(theta,p,K))
          
          
          
          BETA0[[b]][q]<-beta0
          THETA0[[b]][,q]<-theta0
          BETA[[b]][,q]<-beta
         # print(beta1)
        #  print(BETA[[b]][,i])
          THETA[[b]][,,q]<-theta
          
          # main_norm[b]<- (1-alpha)*lambda*sum(norm1)+ alpha*lambda*sum(norm3)
          
          next(b)
          #}
          
        }
        
        
      m_for_beta<-matrix(0,p,max(y))
      for (l in 1:max(y)) {
        m_for_beta[,l]<-BETA[[l]][,q]
      }
      # 
      for (l in 1:nrow(m_for_beta)) {
        med<-mean(m_for_beta[l,])
        m_for_beta[l,]<-m_for_beta[l,]-med
      }
      # 
     # BETA<-lapply(seq_len(max(y)),
      #             function(j)(as.matrix(m_for_beta[,j],nrow = 1,ncol = p)))
      
      
      for (l in 1:max(y)) {
        BETA[[l]][,q]<-m_for_beta[,l]
      }
        
        
        for (x in 1:max(y)) {
          
          beta0=BETA0[[x]][q]
          beta = BETA[[x]][,q]
          # print(beta)
          
          theta =as.matrix(THETA[[x]][,,q]) ; theta_transpose <- t(theta)
          
          
          theta0 = THETA0[[x]][,q]
          
          
          
          
          # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
          n_i[x]<-list(model(beta0, theta0, beta, theta, X, Z))
        }
        E<-matrix(1,N,max(y))
        for (x in 1:max(y)) {
          E[,x]<-exp(unlist(n_i[x]))
          
          
        }
     #   my_pr<-matrix(0,N,(max(y)-1))
        for (x in 1:(max(y))) {
          n_i_b<-matrix(unlist(n_i[x]),N)
          E_1<-E
          E_1[,x]<-0
          pr[x]<-list(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
         # my_pr[,x]<-(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
        }
        #pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
        
        
        #v1_d<-1*(v1==d)
        Dev<-matrix(0,nrow = length(y))
        Dev1<-matrix(0,nrow = length(y))
        for (l in 1:length(Dev)) {
          deviance_y1<-matrix(0,1,max(y))
          deviance_y2<-matrix(0,1,max(y))
          for (d in 1:max(y)) {
            y_d<-1*(y[l]==d)
            # n_i_d<-matrix(unlist(n_i[d]),N)
            prob_d <-matrix(unlist(pr[d]),N)
            # deviance_y1[d]<-y_d*n_i_d[l]
            # deviance_y2[d]<-exp(n_i_d[l])
            deviance_y1[d]<-  Log(y=y_d,pr=prob_d[l]) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
            deviance_y2[d]<-Log(y=y_d,pr=mean(1*(y==d)))
          }
          Dev[l]<-  sum(deviance_y1)
          Dev1[l]<-  sum(deviance_y2)
        }
        # Dev<-matrix(0,nrow = length(y))
        # Dev1<-matrix(0,nrow = length(y))
        # for (l in 1:length(Dev)) {
        #   deviance_y1<-matrix(0,1,max(y))
        #   
        #   deviance_y2<-matrix(0,1,max(y))
        #   for (d in 1:max(y)) {
        #     y_d<-1*(y==d)
        #     # n_i_d<-matrix(unlist(n_i[d]),N)
        #     prob_d <-matrix(unlist(pr[d]),N)
        #     # deviance_y1[d]<-y_d*n_i_d[l]
        #     # deviance_y2[d]<-exp(n_i_d[l])
        #     deviance_y1[d]<-sum(errfun.binomial(y=y_d[l],yhat=prob_d[l],w=rep(1,length(y_d[l])))) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
        #     
        #     deviance_y2[d]=sum(errfun.binomial(y=y_d[l],mean(y_d),w=rep(1,length(y_d[l]))))
        #   }
        #   Dev[l]<-  sum(deviance_y1)
        #   Dev1[l]<-  sum(deviance_y2)
        #   
        #   
        #   #(y,yhat,w=rep(1,length(y)))
        # }
        DEV1[i]<-((2)*sum(Dev1))/length(y)
        DEV[i]<- ((2)*sum(Dev))/length(y) #sum(deviance_y1)# #+sum(main_norm)
        # dev_percent[i]<-(sum(deviance_y2)-sum(deviance_y1))/sum(deviance_y2)
        if(is.null(t)){t=0
        }else{t=t}
        big_t[q]<- t
        Mbeta[q]<-mbeta
        
        Y_hat[q]<- list(matrix(c(y,unlist(pr)),N,(max(y)+1)))
        
        #y_hat<-list(pr)
        
        
        #Y_hat[q]<-y_hat
        
        for (v in 1:max(y)) {
          
          
          b=matrix(0,p)
          
          
          nzero_beta=  BETA[[v]][,q]
          nzero_theta=unlist(THETA[[v]][,,i])
          # for (m in 1:p) {
          #  if (isTRUE(abs(c[m])>0)==T) {b[m]=1
          
          #  }
          
          #}
          non_zero[i,v]<-(sum(nzero_beta!=0 ) )
          my_ntheta<-matrix(0,p)
          # for (w in 1:p ) {
          #   if(isTRUE(sum( nzero_theta[w,])!=0)==T){
          #     my_ntheta[w]<-1 
          #     
          #   }
          # }
          non_zero_theta[i,v]<-(sum( nzero_theta!=0))
          
        }
        
        #non_zero<- data.frame(non_zero,length(b ) )    
       # BETA[q]<-as.list(matrix(unlist(BETA1),p,max(y)))
       # BETA0[q]<-as.list(matrix(unlist(BETA01),1,max(y)))
        #theta01<-theta
       # THETA0[q]<-as.list(matrix(unlist(THETA01),ncol = max(y)))
        
        #sol=list(beta0,beta,theta0,theta)
      #  THETA[q]<-as.list(THETA1)
        Lambda[q]<-lambda
        #Lambda[i]<-max(lambda_i[i,])
       # print(t(active_set1))
        #print(t(active_set2))
        print(c(q,max(non_zero[q,]) ,big_t[q], DEV[q],DEV[q-1]))
        
        
        
        
        
        
        my_q<-my_q+1
        
        # t=new_t/(1+intercept*i)
        
        #t=new_t*exp(1)^(-intercept*i)
        
       # my_v<-my_v+1
        q=q+1
        print(q)
       
        
        
        #if(isTRUE(q==nlambda)==T){
        # q=q
        #   break
        # }else{q=q+1}
      # ########## while loop 
      
      
      
      
    } else if(isTRUE(q>fq)==T){
      
        
     
     i=q
      
    
     lambda<-lambda_i[q]
     
      
      for (b in 1:max(y)) {
        
        lambda<-lambda_i[q]
        
        
        
        beta0<-BETA0[[b]][q-1]
        beta <- BETA[[b]][,q-1]
        # print(beta)
        # print(beta)
        v_beta=BETA[[b]][,q-2]
        theta <-as.matrix(THETA[[b]][,,q-1]) ; theta_transpose <- t(theta)
        v_theta=as.matrix(THETA[[b]][,,q-2]);theta_transpose1 = t(v_theta)
        
        theta0 = THETA0[[b]][,q-1]
        theta01=theta0
        
        beta01=beta0
        beta1 = beta
        theta1 = theta ; theta_transpose1 = t(theta1)
        theta01 = theta
        norm1<-matrix(0,p,1);norm3<-matrix(0,p,1)
        
        
        
        
        
        n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
        
        
        n_i_b<-matrix(unlist(n_i[b]),N)
        E[,b]<-exp(unlist(n_i[b]))
        pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
        pr_b<-matrix(unlist(pr[b]),N)
        
       
        
        
        B<-pr_b*(1-pr_b)
        B[B==0]=as.numeric(10e-9)
        
        y1=1*(y==b)
        
        
        ## Define_Parameters
        #beta01<-beta0
        
        M<-(y1-pr_b)
        
        A<-n_i_b+M/B
        
        
        
        rv<-B*(A-n_i_b)
        r<-(A-n_i_b)
        v=reg(rv,Z)
        
        if (b==max(y)) {
          my_beta0<-matrix(0,(max(y)-1))
          for (l in 1:(max(y)-1)) {
            my_beta0[l]<-BETA0[[l]][i-1]
          }
          beta0<-1-sum(my_beta0)  
          theta0<-matrix(unlist(v[2]))
        } else{
          
          beta0<-matrix(unlist(v[1]))
          theta0<-matrix(unlist(v[2]))
        }
        #r=y1
        
        
        for (v in 1:p) {
          
          # Strong screening rule (https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
          if(i<=1){checkk=2*lambda-lambda_i[1] }else{checkk=( 2*lambda-lambda_i[i-1] )}
          
          
          
          # print(c(i,as.numeric(abs( (t(matrix(X[,j]*c(y) ))%*% matrix((r1) ) )  )/(N*alpha) ), checkk ))
          
          if(as.numeric(abs( (t(matrix( (X[,v]) ) )%*% matrix(((rv))   ) )  )/(N*alpha) )< checkk  ){
            next(v)
          }else if(isTRUE(active_set1[v,b]==1)==T){
            strong_set[v,b]=1
            strong_set1[v,b]=1
            next(v)
          }else {
            
          #  X_j<-X; Z_j<-Z
          #  X_j[,v]<-0
            
           # XZ_theta_j <- lapply(seq_len(p),
           #                      function(j) (matrix(X_j[, j], nrow = N, ncol = K) * Z_j) %*% theta_transpose[, j])
            
            #XZ_theta_j[j]<-0
           # XZ_term_j <- Reduce(f = '+', x = XZ_theta_j)
            
          #  n_j<-as.numeric(beta0)+Z%*%(theta0)+X_j%*%matrix(beta) + XZ_term_j
          #  r_j<-B*(A-n_j)
            
          	
          	#cond1<- as.numeric(abs(t(matrix(X[,v]))%*% matrix(A))/N)
          	
          	
          #	cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[v]))%*%(matrix((A))))/N, alpha*(lambda))), type="F"))
            res_j<- rv+B*model_j(beta[j], theta[j,], X[,j], W, v, Z)
            
          cond1<- as.numeric(abs((t(matrix(X[,v]))%*% matrix( res_j ))/N))
          	
          	
         cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[v]))%*%( res_j ))/N, alpha*(lambda))), type="F"))
          	
          	
          #	cond1<- as.numeric(abs((t(matrix(X[,v]))%*% matrix((y1)))/N))
          	
          	
          #	cond2<-as.numeric(norm(matrix(S_func((t(data.frame(W[v]))%*%(matrix((y1))))/N, alpha*(lambda))), type="F"))
          	
          	
          	#cond1<- as.numeric(abs(t(matrix(B*X[,v]))%*% matrix(r))/N)
          	
          	
          #	cond2<-as.numeric(norm(matrix(S_func((t(B*data.frame(W[v]))%*%(matrix((r))))/N, alpha*(lambda))), type="F"))
            
            #lambda_i[0]<-0
            #if(isTRUE(active_set1[v,b]==1)==T){
            # strong_set[v]<-1  
            # next(v)
            # } else
         if (isTRUE(cond1<= (1-alpha)*(lambda))==T ){
           
           strong_set[v,b]<-0
         } else{
           strong_set[v,b]<-1
         }
         # print(c(cond1,(1-alpha)*(lambda)))
         
         if (isTRUE(cond2<= 2*(1-alpha)*(lambda))==T ){
           
           strong_set1[v,b]<- 0
         } else{
           strong_set1[v,b]<-1
         }
         
            
          }
        }
      
        
      }
      
      
      
      
    while (my_v<=my_V) {
      if(isTRUE(q>mv)==T){
        my_V=ms}else{my_V=my_V}
      print(c(my_v,my_V))
   i=q
    lambda<-lambda_i[q]
    #print(c(q,lambda))
    main_norm<- matrix(0,max(y))
    
    
    
    
    
    
    for (b in 1:max(y)) {
      # lambda_i<- lambda_to_use[,b]
      lambda<-lambda_i[q]
     # strong_set<-matrix(0,p)
     #
      
     # for (it in 1:max_it) {
      
      
      beta0=BETA0[[b]][q-1]
      beta = BETA[[b]][,q-1]
      # print(beta)
      v_beta=beta
      theta =as.matrix(THETA[[b]][,,q-1]) ; theta_transpose <- t(theta)
      v_theta=theta;theta_transpose1 = t(v_theta)
      
      theta0 = THETA0[[b]][,q-1]
      theta01=theta0
      
      beta01=beta0
      beta1 = beta
      theta1 = theta ; theta_transpose1 = t(theta1)
      theta01 = theta
      
      
      
        
      for( iii in 1: max_iter){
      
      n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
      
      
      n_i_b<-matrix(unlist(n_i[b]),N)
      E[,b]<-exp(unlist(n_i[b]))
      pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
      
      # E[,b]<-0
      
      #  pr[b]<-list(matrix(as.numeric(1/( exp(-n_i_b)*rowSums(E)+1),N)))
      pr_b<-matrix(unlist(pr[b]),N)
      
      
      
      
      
      
      B<-pr_b*(1-pr_b)
      
      # O<-B
      # for (m in 1:N) {
      #  if( isTRUE(B[m]==0)==T){
      #    O[m]<- as.numeric(10e-9)
      #  }
      
      # }
      
       B[B==0]=as.numeric(10e-9)
      # B<-O
      
      y1=1*(y==b)
      
      
      ## Define_Parameters
      #beta01<-beta0
      
      M<-(y1-pr_b)
      
      A<-n_i_b+M/B
      
      
      rv<-B*(A-n_i_b)
      iter_prev_score= objective(rv,beta,theta,alpha,lambda)
     # iter_prev_score = objective(rv,beta,theta,alpha,lambda )
      #r<-(A-n_i_b)
      v=reg(rv,Z)
      
      
      if (b==max(y)) {
        my_beta0<-matrix(0,(max(y)-1))
        for (l in 1:(max(y)-1)) {
          my_beta0[l]<-BETA0[[l]][q-1]
        }
        beta0<-1-sum(my_beta0)  
        theta0<-matrix(unlist(v[2]))
      } else{
        
        beta0<-matrix(unlist(v[1]))
        theta0<-matrix(unlist(v[2]))
      }
      #r=y1
      
     
      
      
      DJ_BETA<-matrix(0,1)
      DJ_THETA<-matrix(0,1,K)
      #############################################################
      #############################################################
      #############################################################
      #############################################################
      #############################################################    ########## main loop###########
      
  
        
      
      for (j in 1:p) {
        if(is.null(tt)){
          t=NULL
          
        } else{t=tt[q]}
        
        if(isTRUE(active_set1[j,b]==1)==T & isTRUE(active_set2[j,b]==1)==T){
          
        
          
           #while
          
          #beta1[j]<-beta_new[j];theta1[j,]<-theta_new[j,];t=t
          
          
          value1<-   quadratic(beta,theta, alpha,lambda,beta0,theta0,j,b,W,X,Z,y,N,n_i,xbar,zbar,t=t,r_min_j = res_j)
          #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
          
          beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]])
          
          
          
         
          active_set1[j,b]=1
          active_set2[j,b]=1
          
          next(j)
          
        }else if(as.numeric(abs( (t(matrix(X[,j]))%*% matrix(rv) ) )  )< 2*lambda-lambda_i[(i-1)] ){
          
          next(j)
        }else if(strong_set[j,b]==0 & strong_set1[j,b]==0){
          
         # beta[j]<-0;theta1[j,]<-0
          next(j)
          
        }    else {
          
          
          
          
          # beta1_j<- (N/(t(B)%*%(matrix(X[,j]))^2))*  S_func((t(B*matrix(X[,j]))%*%matrix(r_j))/N,(1-alpha)*lambda)
          
          
          res_j<- rv+B*model_j(beta[j], theta[j,], X[,j], W, j, Z)
          
          # cond3<- norm(matrix(S_func((t(data.frame(W[j]))%*%((r_j)-(B)*matrix(X[,j])*as.numeric(beta1_j))/N),alpha*(lambda))),type="F")
        #  beta_u=beta;beta_new=beta;v_beta=beta;nes=1;n_l=n_i;B_l=B;r_l=r
          #okay_beta=0
         # lam_max<-lambda_i[1]
         # mlambda=lambda_i[i]
        
          beta1_j<- (N/sum(B*(matrix(X[,j])*matrix(X[,j]))))*  S_func((t(matrix(X[,j]))%*%matrix(res_j ))/N,(1-alpha)*lambda)
         #    #  
         #    # 
         #    #       #print(beta1_j)
         #    #  
         #    #       
         #  }
         #  beta1_j<-beta_u[j]
         # # print(beta1_j)
        	
        	
          cond3<- norm(matrix(S_func(t(data.frame(W[j]))%*%( ( res_j )-(B)*matrix(X[,j])*as.numeric(beta1_j))/N,alpha*lambda)),type="F")
         # cond3<- norm(matrix(S_func((t(data.frame(W[j]))%*%(((r+X[,j]*beta[j]+as.matrix(data.frame(W[j]))%*%theta[j,]))-(B)*matrix(X[,j])*as.numeric(beta1_j))/N),alpha*lambda)),type="F")
          
          
          if (isTRUE(cond3<= (1-alpha)*(lambda))==T) {
            beta[j]<-as.numeric(beta1_j)
            #theta[j,]<-as.numeric(0)
            active_set1[j,b]=1
            # v_beta[j]<-unlist(value1[[1]]); v_theta[j,]<- unlist(value1[[2]]); beta[j]<-unlist(value1[[3]]);theta[j,]<-unlist(value1[[4]])
            
            
            next(j)
            
            
          }  else {
            
            
           
            value1<-   quadratic(beta,theta, alpha,lambda,beta0,theta0,j,b,W,X,Z,y,N,n_i,xbar,zbar,t=t,r_min_j = res_j)
            #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
            
            beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]])
            
           
            
           # beta1[j]<-beta_new[j];theta1[j,]<-theta_new[j,];t=t
            
            
            
            
            
            
            active_set1[j,b]=1
            active_set2[j,b]=1
            
            next(j)}
        }
        
        
        
      }
      
      
      
      n_i[b]<-list(model(beta0, theta0, beta, theta, X, Z))
      
      
      n_i_b<-matrix(unlist(n_i[b]),N)
      E[,b]<-exp(unlist(n_i[b]))
      pr[b]<-list(  (matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N)))
      pr_b<-matrix(unlist(pr[b]),N)
      
      
      
      
      
      
      B<-pr_b*(1-pr_b)
      
      
      B[B==0]=as.numeric(10e-9)
      
      
      y1=1*(y==b)
      
      M<-(y1-pr_b)
      
      A<-n_i_b+M/B
      rv<-B*(A-n_i_b)
      tolerance=tol
      if(i>new_t){
        # iter_current_score = objective(
        #   beta0, theta0, beta, theta,
        #   X, Z, y=orig_y,
        #   alpha, lam, W
        # )
        
        iter_current_score= objective(rv,beta,theta,alpha,lambda)
        #print(c(iter_prev_score,iter_current_score,beta))
        if(isTRUE(abs(iter_prev_score - iter_current_score) <tolerance)==T ){
          # print(c(iter_prev_score,iter_current_score))
          break  # Converged on lam_i
          
        }else{next(iii)}
        
      }else{break}
      
      
      }#maxt it
      
    
        BETA01[b]<-list(beta0) 
        
        BETA1[b]<-list(beta)
        
        THETA01[b]<-list(matrix(theta0,1,K))
        THETA1[b]<-list(as.matrix(theta,p,K))
        
        
        
        BETA0[[b]][q]=beta0
        THETA0[[b]][,q]=theta0
        BETA[[b]][,q]=beta
       # print(beta1)
       # print(BETA[[b]][,i])
        THETA[[b]][,,q]=theta
        
       # main_norm[b]<- (1-alpha)*lambda*sum(norm1)+ alpha*lambda*sum(norm3)
        
        next(b)
        #}
        
      }
    
    
    m_for_beta<-matrix(0,p,max(y))
    for (l in 1:max(y)) {
      m_for_beta[,l]<-BETA[[l]][,q]
    }
    # 
    for (l in 1:nrow(m_for_beta)) {
      med<-mean(m_for_beta[l,])
      m_for_beta[l,]<-m_for_beta[l,]-med
    }
    # 
   # BETA<-lapply(seq_len(max(y)),
    #             function(j)(as.matrix(m_for_beta[,j],nrow = 1,ncol = p)))
    
    
    for (l in 1:max(y)) {
      BETA[[l]][,q]<-m_for_beta[,l]
    }
    
    
    for (x in 1:max(y)) {
      
      beta0=BETA0[[x]][q]
      beta = BETA[[x]][,q]
      # print(beta)
      
      theta =as.matrix(THETA[[x]][,,q]) ; theta_transpose <- t(theta)
      
      
      theta0 = THETA0[[x]][,q]
      
      
      
      
      
      
      # n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
      n_i[x]<-list(model(beta0, theta0, beta, theta, X, Z))
    }
    E<-matrix(1,N,max(y))
    for (x in 1:max(y)) {
      E[,x]<-exp(unlist(n_i[x]))
      
      
    }
    my_pr<-matrix(0,N,(max(y)-1))
    for (x in 1:(max(y))) {
      n_i_b<-matrix(unlist(n_i[x]),N)
      E_1<-E
      E_1[,x]<-0
      pr[x]<-list(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
    #  my_pr[,x]<-(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
    }
   # pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
    
    #v1_d<-1*(v1==d)
    Dev<-matrix(0,nrow = length(y))
    Dev1<-matrix(0,nrow = length(y))
    for (l in 1:length(Dev)) {
      deviance_y1<-matrix(0,1,max(y))
      deviance_y2<-matrix(0,1,max(y))
      for (d in 1:max(y)) {
        y_d<-1*(y[l]==d)
        # n_i_d<-matrix(unlist(n_i[d]),N)
        prob_d <-matrix(unlist(pr[d]),N)
        # deviance_y1[d]<-y_d*n_i_d[l]
        # deviance_y2[d]<-exp(n_i_d[l])
        deviance_y1[d]<-  Log(y=y_d,pr=prob_d[l]) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
        deviance_y2[d]<-Log(y=y_d,pr=mean(1*(y==d)))
      }
      Dev[l]<-  sum(deviance_y1)
      Dev1[l]<-  sum(deviance_y2)
    }
    # Dev<-matrix(0,nrow = length(y))
    # Dev1<-matrix(0,nrow = length(y))
    # for (l in 1:length(Dev)) {
    #   deviance_y1<-matrix(0,1,max(y))
    #   
    #   deviance_y2<-matrix(0,1,max(y))
    #   for (d in 1:max(y)) {
    #     y_d<-1*(y==d)
    #     # n_i_d<-matrix(unlist(n_i[d]),N)
    #     prob_d <-matrix(unlist(pr[d]),N)
    #     # deviance_y1[d]<-y_d*n_i_d[l]
    #     # deviance_y2[d]<-exp(n_i_d[l])
    #     deviance_y1[d]<-sum(errfun.binomial(y=y_d[l],yhat=prob_d[l],w=rep(1,length(y_d[l])))) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
    #     
    #     deviance_y2[d]=sum(errfun.binomial(y=y_d[l],mean(y_d),w=rep(1,length(y_d[l]))))
    #   }
    #   Dev[l]<-  sum(deviance_y1)
    #   Dev1[l]<-  sum(deviance_y2)
    #   
    #   
    #   #(y,yhat,w=rep(1,length(y)))
    # }
    DEV1[i]<-((2)*sum(Dev1))/length(y)
    DEV[i]<-((2)*sum(Dev))/length(y) #sum(deviance_y1)# #+sum(main_norm)
    # dev_percent[i]<-(sum(deviance_y2)-sum(deviance_y1))/sum(deviance_y2)
    BETA01_past<-BETA01
    
    BETA1_past<-BETA1
    
    THETA01_past<-THETA01
    THETA1_past<-THETA1
    if(is.null(t)){
      t=0
    }else{t=t}
     big_t[q]<- t 
    Mbeta[q]<-mbeta
    
    Y_hat[q]<- list(matrix(c(y,unlist(pr)),N,(max(y)+1)))
    
    #y_hat<-list(pr)
    
    
    #Y_hat[q]<-y_hat
    
    for (v in 1:max(y)) {
      
      
      b=matrix(0,p)
      
      
      nzero_beta=  BETA[[v]][,q]
      nzero_theta=unlist(THETA[[v]][,,i])
      # for (m in 1:p) {
      #  if (isTRUE(abs(c[m])>0)==T) {b[m]=1
      
      #  }
      
      #}
      non_zero[i,v]<-(sum(nzero_beta!=0 ) )
      my_ntheta<-matrix(0,p)
      # for (w in 1:p ) {
      #   if(isTRUE(sum( nzero_theta[w,])!=0)==T){
      #     my_ntheta[w]<-1 
      #     
      #   }
      # }
      non_zero_theta[i,v]<-(sum( nzero_theta!=0))
      
    }
    
    #non_zero<- data.frame(non_zero,length(b ) )    
   # BETA[q]<-as.list(matrix(unlist(BETA1),p,max(y)))
  #  BETA0[q]<-as.list(matrix(unlist(BETA01),1,max(y)))
    #theta01<-theta
    #THETA0[q]<-as.list(matrix(unlist(THETA01),ncol = max(y)))
    
    #sol=list(beta0,beta,theta0,theta)
   # THETA[q]<-as.list(THETA1)
    Lambda[q]<-lambda
    #Lambda[i]<-max(lambda_i[i,])
    print(c(q,max(non_zero[q,]) ,big_t[q], DEV[q],DEV[q-1]))
    
    
    
    
    
    
    my_q<-my_q+1
    
    # t=new_t/(1+intercept*i)
    
    #t=new_t*exp(1)^(-intercept*i)
    
    my_v<-my_v+1
    q=q+1
    print(q)
    if( isTRUE(q>nlambda)==T){
      #my_ok<-1
      print(q)
      break
    }
    } 
    } ## i>1 
     
   
   
    
   
    
    my_v<-sv
    my_V<-for_v
    q=q
    
    sumact<-matrix(0,max(y))
    for (nn in 1:max(y)) {
      act<-(active_set1[,nn])
      sact<-(sec_active1[,nn])
      if(isTRUE(all(act==sact)==T))
        sumact[nn]<-1
    }
    
    
  tolerance=tol
      
   if(isTRUE(q>st)==T){
    if(isTRUE(max(sumact)>=1)==T | isTRUE(q>nlambda)==T){
      my_ok<-1
      break
      
    }else if ( isTRUE(sum(sumact)<1)==T & isTRUE(q<=nlambda)==T){
      sec_active1<-active_set1
      sec_active2<-active_set2
      my_ok<-0
      #next(i)
    }else{
      my_ok<-0
    }
   } else if(isTRUE(q>nlambda)){
     my_ok=1}else{my_ok=0}
    #if(all(colSums(active_set1)==colSums(sec_active1))==T){
    #  my_ok<-1
      
   # }
    
 # } #### i
  } ## for ok
  
  
  
  
  dev_percent<-matrix(0,nlambda)
  for (i in 1:nlambda) {
    # if (i==1) {dev_percent[i]<-0
    
    #} else {
    dev_percent[i]<-1-((DEV[i])/DEV1[i])
    
    # }
    # dev_percent[i]<-((DEV[i])/DEV[i-1])
    
  }
  #}
  nzero<-matrix(0,nrow(ACTIVE))
  nzero_int<-matrix(0,nrow(ACTIVE1))
  for (s in 1:nlambda) {
    # nzero[s]<-max(ACTIVE[s,]) 
    nzero_int[s]<-max(non_zero_theta[s,])
    nzero[s]<-max(non_zero[s,])
    
  }
  
  
  pred<-data.frame(Lambda=matrix(Lambda,ncol = 1),nzero=nzero,nzero_int=nzero_int,DEV=matrix(DEV,ncol=1),nullDEV=matrix(DEV1,ncol=1),Dev_rat=dev_percent)
  
  
  return(list(beta0=BETA0,beta=BETA,theta0=THETA0,theta=THETA,y_hat=Y_hat,path=pred,Lambdas=Lambda,big_T=big_t,Mbeta=Mbeta,non_zero=non_zero))
}



  
  
  

  
  
  
  
  
  

############################################
############################################
############################################
############################################
##################Cross validation#########


predict_lasso<-function(object ,X,Z,y,lambda=NULL){
  lambda.arg=lambda
  if(is.null(lambda.arg))
  { lambda=object$Lambdas;isel=1:length(lambda)}
  
  if(!is.null(lambda.arg)){
    
    isel=as.numeric(knn1(matrix(object$Lambdas,ncol=1),matrix(lambda.arg,ncol=1),1:length(object$Lambdas)))
    
  }
  
  
  #print(isel)
  
  N <- nrow(X)
  
  p <- ncol(X)
  
  K <- ncol(Z)
  
  
  yh=array(0,length(isel))
  DEV=matrix(NA,length(isel))
  my_theta<-array(0,c(ncol(X),ncol(Z), length(isel)))
  
  
 #pBETA0<-matrix(0,nrow = length(isel)); pBETA<-matrix(0,nrow = length(isel)); pTHETA0<-matrix(0,nrow = length(isel)); pTHETA<-matrix(0,nrow =length(isel))
  
  
  pBETA0<-lapply(seq_len(max(y)),
                function(j)(matrix(0,nrow = length(isel))))
  
  pBETA<-lapply(seq_len(max(y)),
                function(j)(matrix(0,nrow=p,ncol=length(isel))))
   
  pTHETA0<-lapply(seq_len(max(y)),
                 function(j)(matrix(0,nrow=K,ncol=length(isel))))
  
  pTHETA<-lapply(seq_len(max(y)),
                function(j)(array(0,c(p,K,length(isel)))))
  
  
  
  iii=0
  for(m in isel){
    iii=iii+1
    
   # pred_beta<-lapply(seq_len(max(y)),
    #                  function(j)(matrix(0,nrow = 1,ncol = ncol(X))))
    #pred_theta<-lapply(seq_len(max(y)),
   #                    function(j)(matrix(0,nrow = p,ncol = K)))
  #  pred_theta0<-lapply(seq_len(max(y)),
    #                    function(j)(matrix(0,nrow = 1,ncol = K)))
    
  #  pred_beta0<-lapply(seq_len(max(y)),
     #                  function(j)(0))
    
    
    z=m
    #BETA0<-matrix(unlist(object$beta0[z]),1,max(y))
    
    
   # BETA<-as.data.frame(object$beta[z])
   # THETA0<-as.data.frame(object$theta0[z])
   # THETA<-object$theta[z]
    
    
    
    
    n_i<-lapply(seq_len(max(y)),
                function(j)(matrix(0,nrow = N)))
    pr<-lapply(seq_len(max(y)),
               function(j)(matrix(0,nrow = N)))
    
    
    for (x in 1:max(y)) {
      #theta<- matrix(unlist(THETA[[1]][x]),p,K)
      
      
      beta0<-object$beta0[[x]][z]
      beta <- object$beta[[x]][,z]
    
      theta <- as.matrix(object$theta[[x]][,,z]) ; theta_transpose <- t(theta)
      
      
      theta0 <- object$theta0[[x]][,z]
      
      
      pBETA0[[x]][iii] <-beta0
      pBETA[[x]][,iii] <-beta
      pTHETA[[x]][,,iii] <-theta
      pTHETA0[[x]][,iii] <-theta0
      
      
      
      
     
      
     # pred_theta[x]<-list(as.matrix(theta,p,K))
     # pred_beta0[x]<-beta0
     # pred_theta0[x]<-theta0
     # pred_beta[x]<-beta
     # beta=matrix(unlist(BETA[x]),p)
      #n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
      n_i[x]<-list(model(beta0, theta0, beta, theta, X, Z))
      
    }
    E<-matrix(0,N,max(y))
    for (x in 1:max(y)) {
      E[,x]<-exp(unlist(n_i[x]))
    }
    for (x in 1:max(y)) {
      n_i_b<-matrix(unlist(n_i[x]),N)
      
      pr[x]<-list(matrix(as.numeric(exp(n_i_b)/( (rowSums(E)))),N))
      
    }
    
    Dev<-matrix(0,nrow = length(y))
    
    for (l in 1:length(Dev)) {
      deviance_y1<-matrix(0,1,max(y))
      deviance_y2<-matrix(0,1,max(y))
      for (d in 1:max(y)) {
        y_d<-1*(y[l]==d)
        # n_i_d<-matrix(unlist(n_i[d]),N)
        prob_d <-matrix(unlist(pr[d]),N)
        # deviance_y1[d]<-y_d*n_i_d[l]
        # deviance_y2[d]<-exp(n_i_d[l])
        deviance_y1[d]<-  Log(y=y_d,pr=prob_d[l]) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
      }
      Dev[l]<-  sum(deviance_y1)
      
    }
    
    # Dev<-matrix(0,nrow = length(y))
    # 
    # for (l in 1:length(Dev)) {
    #   deviance_y1<-matrix(0,1,max(y))
    #   
    #   deviance_y2<-matrix(0,1,max(y))
    #   for (d in 1:max(y)) {
    #     y_d<-1*(y[l]==d)
    #     # n_i_d<-matrix(unlist(n_i[d]),N)
    #     prob_d <-matrix(unlist(pr[d]),N)
    #     # deviance_y1[d]<-y_d*n_i_d[l]
    #     # deviance_y2[d]<-exp(n_i_d[l])
    #     deviance_y1[d]<-sum(errfun.binomial(y=y_d,yhat=prob_d[l],w=rep(1,length(y_d)))) #+(1-y_d)*Log(y=(1-y_d),pr=(1-prob_d[l]))
    #     
    #    # deviance_y2[d]=sum(errfun.binomial(y=y_d[l],mean(y_d),w=rep(1,length(y_d[l]))))
    #   }
    #   Dev[l]<-  sum(deviance_y1)
    #   #Dev1[l]<-  sum(deviance_y2)
    #   
    #   
    #   #(y,yhat,w=rep(1,length(y)))
    # }
    #DEV1[i]<-((2)*sum(Dev1))/length(y)
    DEV[iii]<-((2)*sum(Dev))/length(y)
    
   # DEV[ii]<- ((2)*sum(Dev))#+
    yh[iii]<- list(matrix(unlist(pr),N,max(y)))
   
    
    
  }
  
  
  
  return(list(y_hat=yh,beta0=pBETA0,beta=pBETA,theta0=pTHETA0,theta=pTHETA,deviance=DEV))
}



