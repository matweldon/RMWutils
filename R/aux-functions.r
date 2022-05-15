#' A set-seed clone that initialises at a random value, and returns the value
#'
#' @param x an optional seed to provide
#'
#' @return the value of the seed that has been set
#'
#' @export
mysetseed <- function(x=NULL){
  if(exists("given.seed")){
    myseed <- given.seed
  }else if(!is.null(x)){
    myseed <- x
  }else{
    #Time seed
    mytime <- as.numeric(Sys.time())
    myseed <- round((mytime-floor(mytime))*1e5)
  }
  while(myseed<1e4){ ### seed is always five digits (for filename formatting)
    myseed <- myseed*10
  }
  set.seed(myseed)
  return(myseed)
}

#' Get p-spline trajectory with conf bands
#'
#' @export

yhat <- function(X,vc,pars,z=c(-1.96,0,1.96),FUN=NULL){
  if(is.null(FUN)){FUN <- function(x){x}}
  #getf <- get(FUN)
  se <- sqrt(diag(X %*% vc %*% t(X)))
  lp <- X %*% pars
  lpm <- cbind(lp,lp,lp) + (as.matrix(se) %*% z)
  ym <- FUN(lpm)

  attr(ym,"lpm") <- lpm
  attr(ym,"se") <- se
  return(ym)
}

#' Get Hessian for optimisation without errors
#'
#' @export
robustHess <- function(par,FUN,...,sd=T){
  tryCatch(
    {hss <- optimHess(par,FUN,...) ;
    if(sd==T){
      res <- sqrt(diag(solve(-hss)));
      if(length(res!=length(par))){return(rep(NA,length(par)))
      }else{return(res)}
    }
    else{return(solve(-hss))}
    }
    #,warning= function(w){}
    ,error= function(e){warning("Errors produced");
      if(sd==T){return(rep(NA,length(par)))
      }
      else{
        return(matrix(NA,length(par),length(par)))
      }
    }
  )
}




#' Numerically safe log(sum(exp(x)))
#'
#' Implements the log-sum-exp trick
#'
#' @param x a numeric vector
#'
#' @export
logsumexp <- function(x){
  max.x <- max(x)
  x2 <- x-max.x
  lse <- log(sum(exp(x2)))+max.x
  return(lse)
}


#' Convert vector of indices to logical and vice versa
#'
#' @export
idx2tf <- function(x,n=max(x)){1:n %in% x}
#' @export
tf2idx <- function(x){(1:length(x))[x]}

#' Insert into a vector
#'
#' @export
insert <- function(v,v1,after=length(v)){
  ### Identical to append except allows vector "after"
  ### After can be a scalar or vector
  ### If it is a scalar all the elements of v1 are inserted at the same point
  v0 <- c(v,v1)
  if(length(v1)>1&length(after)==1){after <- rep(after,length(v1))}
  ids <- c(1:length(v),after+0.1)
  return(v0[order(ids)])
}



#' Create a matrix of dummy variables
#'
#' @export
make_dummy <- function(n,reps,ref=F){

  if(ref==T){
    dum <- list(diag(n)[,-1])
  }else{
    dum <- list(diag(n))
  }
  return(do.call(rbind,rep(dum,reps)))
}




##----------------------------------------------------------------
##                     Effective dimension of spline
##----------------------------------------------------------------

#' The effective degrees of freedom of a fitted spline
#'
#' The effective number of parameters of a fitted spline depends
#' on the degree of the spline, the number of knots, and the strength
#' and degree of the prior used to regularise the spline.
#'
#' @export
eff.dim <- function(sd=NULL,B,lam=NULL,d=2,k=NULL,scale=1){
  if(is.null(k)){k <- ncol(B)}
  if(length(scale)==1){scale <- rep(scale,k-1)}
  if(is.null(lam)){
    if(is.null(sd)){
        stop("Either sd or lam must be supplied")
      }else{
        lam <- 1/(2*sd^2)
    }
  }
  if(d>1){
    D <- diff(diag(k-1),diff=d-1)%*%diag(scale)%*%diff(diag(k),diff=1)
  }else{D <- diag(scale)%*%diff(diag(k),diff=1)}

  H <- crossprod(B)%*%solve(crossprod(B) + lam*crossprod(D))

  return(sum(diag(H)))
}

#' Effective degrees of freedom of a matrix with prior
edf <- function(xtx,Q){
  H <- xtx%*%solve(xtx + Q)
  return(sum(diag(H)))
}

# Can't remember what this is?
ed.target <- function(ed,B,d=2,k=NULL,scale=1){
  if(is.null(k)){k <- ncol(B)}

  obj <- function(sd){
    ed1 <- eff.dim(sd=sd,B=B,d=d,k=k,scale=scale)
    return((ed1-ed)^2)
  }
  sdobj <- optimise(obj,interval=c(0.001,100))
  return(sdobj$minimum)
}


###---------------------------------------------------------------
###         Make a legend outside multiple plots
###---------------------------------------------------------------
### Remember to set par(oma=c(...)) with large margin where you want the plot
outer_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

