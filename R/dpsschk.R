#'  function to calculate tapers
#'
#' Called by spectral estimation routines to generate the tapers
#'
#' @param N (number of samples
#' @param K - number of tapers, often 2*nw for  spectrum estimation purposes.
#' @param NW - A positive double-precision number, the time-bandwidth parameter
#' @param Fs sampling frequency
#' @return tapers   
#' 
#' @export
dpsschk <-function(N,k,nw, Fs){
    tapers <-multitaper::dpss(N=N,k=k,nw=nw)
    tapers$v = tapers$v*sqrt(Fs)
    tapers$v
}
