#' Multi-taper coherency - point process times
#'
#' Multi-taper coherency - point process times
#'
#' @param  data1  array of spike times 
#' @param  data2  array of spike times
#' @param Fs sampling frequency
#' @param fpass (fmin, fmax), frequency band to be used in the calculation. Defaults to (0, Fs/2)
#' @param pad   (padding factor for the FFT). Not currently used, i.e. no padding is done
#' @param nw   A positive double-precision number, the time-bandwidth paramter for the tapers
#' @param k   A positive integer, the number of tapers, defaults to  2*nw-1
#' @param fscorr  finite size corrections -- not implemented yet
#'
#' @return  S1     the estimated spectrum for data1
#' @return  S2     the estimated spectrum for data2
#' @return  C12    the coherency
#' @return  f      the frequencies
#' 
#' @export

coherencypt<-function(data1,data2,Fs  =   3000,fpass  =  c(0, 1500),  pad  =  0,nw=50 ,k= 2*nw-1){
    mintime1 <- min(data1)
    maxtime1 <- max(data1)
    mintime2 <- min(data1)
    maxtime2 <- max(data1)
    
    mintime <-min(mintime1,mintime2)
    maxtime<-max(maxtime1,maxtime2)
    dt=1/Fs;
    t <- seq(from=mintime-dt,by=dt,to=maxtime+dt) # time grid for prolates
    N<-length(t)# % number of points in grid for dpss
    nfft <- max(2^(pracma::nextpow2(N)+pad),N) #number of points in fft of prolates
    aa<-getfgrid(Fs,nfft,fpass)
    f<-aa$f
    findx<-aa$findx

    tapers <- dpsschk(N,k,nw, Fs) #check tapers
    temp1 <-mtfftpt(data1,tapers,nfft,t,f,findx)  #; % mt fft for point process times
    temp2 <-mtfftpt(data2,tapers,nfft,t,f,findx)  #; % mt fft for point process times
    J1<-temp1$J 
    J2<-temp2$J 
    S1<-apply(Conj(J1)*J1,1,mean)
    S2<-apply(Conj(J2)*J2,1,mean)
    S12<-apply(Conj(J1)*J2,1,mean)

    C12 <- abs(S12/sqrt(S1*S2))
    
    aa<-list()
    aa$S1 <- S1
    aa$S2 <- S2
    aa$C12 <- C12
    aa$f  <- f
    aa
}
