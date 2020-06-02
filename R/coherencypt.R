#' Multi-taper coherency - point process times
#'
#' Multi-taper coherency - point process times
#'
#' @param  data1  array of spike times 
#' @param  data2  array of spike times 
#' @param  params
#' @param fscorr  finite size corrections
#' @param  t   time grid over which the tapers are to be calculated
#'
#'  @return  C          frequencies)
#'  @return findx     index of the frequencies in the full frequency grid
#' 
#' @export
coherencypt<-function(data1,data2,Fs  =   3000,err  =  c(2, 0.0500),fpass  =  c(0, 1500),
                      pad  =  0,tt =  matrix(c(50, 99),1,2),tapers=matrix(c(50, 99),1,2),fscorr =  0){
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

    tapers <- dpsschk(tapers,N,Fs) #check tapers
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
