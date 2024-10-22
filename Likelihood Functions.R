##############################################################################
#  Title: Tuberculosis Transmission Dynamics in a High Incidence Setting     #
###############################################################################
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' Functions to estimate R and k
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' _______________________________________________________________________________________________
##' Likelihood Function
##' For use in parameter estimation 
##'      @param Y 3-column data frame or matrix containing 
##'                  [1] Custer Size
##'                  [2] Index Cases
##'                  [3] Censored status
##'      @param R NB R value to be plugged into the likelihood
##'      @param k NB k value to be plugged into the likelihood
##'      - - - - - - - - - - - - -       
##'      @return Sum of the log-likelihoods
##' _______________________________________________________________________________________________

likelihood <- function(Y,R,k) {
  p_function <- function(y,n){         
    exp(log(n)-log(y)+lgamma(k*y+y-n)-(lgamma(k*y)+lgamma(y-n+1))+(y-n)*log(R/k)-(k*y+y-n)*log(1+R/k))
  }
  ya <- Y[Y[,3]==0,] 
  yb <- Y[Y[,3]==1,]
  
  liks_a <- log(p_function(ya[,1],ya[,2])) 
  
  liks_b <- numeric(nrow(yb))              
  if(nrow(yb) > 0){                        
    for (i in 1:nrow(yb)){
      y <- yb[i,1]
      n <- yb[i,2]
      if (y==1){                           
        liks_b[i] <- 0                     
      } else{
        liks_b[i] <- log(max(10^-300, 1 - sum(p_function(1:(y-1),n)), na.rm = TRUE))  
      }}}
  sumliks <- sum(liks_a,liks_b)
  return(sumliks)
  #return(cbind(liks_a,liks_b))
}

##' _______________________________________________________________________________________________
##' Surface likelihood function
##' Calculates likelihoods over a range of R and k values
##'      @param data 3-column data frame or matrix containing 
##'                  [1] Custer Size
##'                  [2] Index Cases
##'                  [3] Censored status
##'      @param Rrange Range of R values
##'      @param krange Range of k values
##'      - - - - - - - - - - - - -       
##'      @return Rrange by krange Matrix with likelihoods
##' _______________________________________________________________________________________________

surflike <- function(data, Rrange, krange){
  likesurf <- matrix(NA, nrow = length(Rrange),length(krange))
  for(i in 1:length(Rrange)){
    for(j in 1:length(krange)){
      likesurf[i,j] <- likelihood(data,Rrange[i],krange[j])
    }
  }
  return(likesurf)
}

##' _______________________________________________________________________________________________
##' Parameter Estimation
##' Estimates MLE and confidence interval for R and k
##'      @param ls likelihood surface data
##'      @param ls_max logical likelihood surface data identifying max  (ls_max <- ls==max(ls))
##'      @param conf.interval Desired confidence interval (as decimal, i.e. 0.95)
##'      - - - - - - - - - - - - -       
##'      @return Point, lower, and upper bound estimates for R and k 
##' _______________________________________________________________________________________________

calc_profile <- function(ls, ls_max, conf.interval = 95, Rrange = Rrange, krange = krange){
  chiV <- qchisq(conf.interval/100, df = 1)/2
  prfk <- apply(ls,2,function(x){max(x)})
  prfk2 <- krange[prfk-max(prfk)>-chiV]
  prfR <- apply(ls,1,function(x){max(x)})
  prfR2 <- Rrange[prfR-max(prfR)>-chiV]
  
  output <- rbind(cbind(Rrange[sum(seq(1,length(Rrange))%*%ls_max)],min(prfR2),max(prfR2)),
                  cbind(krange[sum(ls_max%*%seq(1,length(krange)))],min(prfk2),max(prfk2)))
  colnames(output) <- c("point_est","lower_ci","upper_ci")
  rownames(output) <- c("R","k")
  return(output)
}

