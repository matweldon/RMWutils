#' Construct a posterior function for optimisation
#'
#' Returns a function of a single parameter vector
#'
#' @param FUNS list of functions (log-likelihoods and log-priors)
#' @param master logical matrix giving indices of par vector for each function
#' @param args List of arguments to each function
#'
#' @return a function that accepts a single parameter vector and
#'         evaluates the log prior. The data supplied in 'args' is
#'         "baked in" to the function environment and does not need
#'         to be resupplied.
#'
#' @export
postFUN <- function(FUNS,master,args){

  if(length(FUNS)!=ncol(master) | ncol(master)!=length(args)){
    stop("Some input to postFUN is the wrong size.")
  }

  FUNSS <- list()
  for(i in 1:length(FUNS)){
    fun1 <- function(){
      force(args1 <- args[[i]])
      force(idx1 <- master[,i])
      force(ff <- FUNS[[i]])
      return(function(x){do.call(ff,c(list(par=x[idx1]),args1))})
    }
    FUNSS[[i]] <- fun1()
  }


  plus <- function(f1,f2){force(f1);force(f2);return(function(x){f1(x)+f2(x)})}
  post <- Reduce(plus,FUNSS)
  return(post)
}


##----------------------------------------------------------------
##                     Priors
##----------------------------------------------------------------

#' A normal prior with a mix of L1 and L2 penalties
#'
#' @export
norm_prior <- function(par,sdf=3,l1penal=0,idx=1:length(par)){

  l2penal <- 1/(2*sdf^2)

  theta <- par[idx]

  ridge <- -l2penal*sum(theta^2)
  lasso <- -l1penal*sum(abs(theta))

  return(ridge+lasso)
}


#' A gamma prior
#'
#' @export
gamma_prior <- function(par,shape=1.5,scale=1
                        ,idx=1:length(par)
                        ,trans.exp=F
){
  theta <- par[idx]
  sig <- abs(theta)*(1-trans.exp)+exp(theta)*trans.exp

  ### I'm not at all sure this needs a jacobean, since both
  ### the prior and the likelihood are on the same scale.
  jacob <- 0*(1-trans.exp)+exp(theta)*trans.exp
  return(sum(dgamma(sig,shape=shape,scale=scale,log=T))+jacob)
}

#' A prior for school means
#'
#' @export
sch_mu_prior <- function(par,sdf=10,idx=1:length(par)){
  sum(dnorm(par[idx],sd=sdf,log=T))
}


#' A prior on d'th differences for p-spline models
#'
#' @export
spline_prior <- function(par,d=2,lambda=1,zero=T
                         ,scale=1,idx=1:length(par)){
  ### P-spline prior
  ### Scale is vector of length k-d for scaling differences

  par1 <- par[idx]
  k <- length(par1)#+ zero

  #B <- par1[idx_B]
  if(zero){B<- c(0,par1);k<-k+1 }else{B <- par1}

  #   if(d>=k){
  #     stop("d must be < number of B-spline co-efficients")
  #   }

  ### Matrix of d-differences
  D <- diff(diag(k),differences=d)

  return(-lambda*sum(((D%*%B)/scale)^2))
}

#' Another version of a p-spline prior
#'
#' @export
spline2_prior <- function(par,sigma=1 #,zero=T
                          ,scale=rep(1,length(par)),idx=1:length(par)){
  ### P-spline prior
  ### Only works for d=2
  ### Scale is vector of length k-1 (+1 for intercept) for scaling differences

  par1 <- par[idx]
  k <- length(par1)#+ zero

  #B <- par1[idx_B]
  #if(zero){B<- c(0,par1);k<-k+1 }else{B <- par1}
  B <- c(0,par1)
  k <- k+1

  J <- diff(diag(k-1),diff=1)%*%diag(scale)%*%diff(diag(k),diff=1)


  lambda <- 1/(2*sigma^2)
  return(-lambda*sum((J%*%B)^2))
}
