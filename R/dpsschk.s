dpsschk <-function(tapers,N,Fs){
#'  function to calculate tapers
#'
#' Called by spectral estimation routines to generate the tapers
#'
#' @param tapers  NW K - time-bandwidth product, number of tapers
#' @param N (number of samples
#' @param Fs sampling frequency
#' @return tapers   
#' 
#' @export

    sz<- dim(tapers)
    if (sz[1]==1 && sz[2]==2){
        tapers <-dpss(N,tt[2],tt[1])
        tapers$v = tapers$v*sqrt(Fs)
    }
    tapers$v
}
