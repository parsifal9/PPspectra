getfgrid <-function(Fs,nfft,fpass){
#'  function that gets the frequency grid associated with a given fft based computation
#'
#' Called by spectral estimation routines to generate the frequency axes
#'
#' @param Fs sampling frequency associated with the data
#' @param nfft number of points in fft
#' @param fpass band of frequencies at which the fft is being calculated (fmin fmax) in Hz
#' @return f         frequencies)
#' @return findx     index of the frequencies in the full frequency grid
#' 
#' @export

    a<-list()
    df <- Fs/nfft
    f <- seq(from=0,to=Fs,by=df)# all possible frequencies
    f <- f[1:nfft]
    findx <- which(f>=fpass[1] & f<=fpass[length(fpass)])
    f <- f[findx]
    a$f <- f
    a$findx <-findx
    a
}
