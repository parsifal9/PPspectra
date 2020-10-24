#'  utility function
#'
#' returns  (exp(-1i* aa %*% t(bb)))%*% cc, The same as  exp(-1i*outer(aa,bb))%*% cc but the outer product aa %*% t(bb) may be very large
#'
#' @param aa vector
#' @param bb vector
#' @param cc matrix
#'
#' @return J matrix   
#' 
f1<-function(aa,bb,cc){
    qq<-length(bb)
    J<-matrix(0,length(aa),dim(cc)[[2]])
    for ( mm in 1:length(aa)){
        temp1<-aa[mm]*bb
        exponential <- exp(-1i*temp1)
        J1 <- exponential %*%cc
        J[mm ,] <-J1
    }
    J
}
