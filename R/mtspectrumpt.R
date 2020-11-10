#'  multi-taperd spectrum of data preseted as spike locations
#'
#' multi-taperd spectrum of data preseted as spike locations
#'
#' @param PP the data
#' @param Fs sampling frequency
#' @param fpass (fmin, fmax), frequency band to be used in the calculation. Defaults to (0, Fs/2)
#' @param pad   (padding factor for the FFT). Not currently used, i.e. no padding is done
#' @param nw   A positive double-precision number, the time-bandwidth paramter for the tapers
#' @param k   A positive integer, the number of tapers, defaults to  2*nw-1
#'
#' @return S  the estimated spectrum
#' @return f  the frequencies
#' 
#' @export
mtspectrumpt<-function(PP, Fs  =   3000,fpass  =  c(0, Fs/2),pad  =  0,nw=50 ,k= 2*nw-1,fscorr =  0){
    mintime <- min(PP)
    maxtime <- max(PP)
    dt <- 1/Fs #% sampling time
    t <- seq(from=mintime-dt,by=dt,to=maxtime+dt)
    
    N<-length(t)# % number of points in grid for dpss
    nfft <- max(2^(pracma::nextpow2(N)+pad),N) #number of points in fft of prolates
    
    aa<-getfgrid(Fs,nfft,fpass)
    f<-aa$f
    findx<-aa$findx
    
    tapers <- dpsschk(N,k,nw, Fs)
    
    temp <-mtfftpt(PP,tapers,nfft,t,f,findx)  #; % mt fft for point process times
    
    ## C <-1
    ## K <-dim(tapers)[2]
    ## nfreq <- length(f); #% number of frequencies
    
    ## #H<-apply(tapers$v,2,fft)
    ## temp<-matrix(0, 1269,99)
    ## H<-apply(rbind(tapers,temp),2,fft)
    
    
    ## # dim(A) 10 20  -- apply(A,1,sum) means sum the rows. apply(A,2,sum) means sum the columns
    ## # However in matlab  fft(A,n,1) means apply fft along the first dimension i.e. to the columns. 
    
    
    ## #H <- fft(tapers$v) #% fft of tapers  fft(X,n,1); length of X is less than n, then X is padded with trailing zeros to length n.
    ## dim(H) #[1] ] 2827   99
    ## H <- H[findx,] #; % restrict fft of tapers to required frequencies
    ## # 2049   99
    
    ## #H=fft(tapers,nfft,1);  % fft of tapers
    ## w <- 2*pi*f    #; % angular frequencies at which ft is to be evaluated
    ## sp <- 0
    ## Msp<- 0
    ## Nsp<- 0
    ## ch <- 1
    ## dtmp <- S1
    
    ## indx <- which(dtmp>=min(t) & dtmp<=max(t))
    ## dtmp <- dtmp[indx]
    ## Nsp[ch] <- length(dtmp)
    ## Msp[ch]  <- Nsp[ch]/length(t);
    
    ## data_proj <-apply(tapers,2,fm<-function(x){pracma::interp1(t,x,dtmp) })
    
    ## aa<-outer(w,(dtmp-t[1]))
    ## exponential <- exp(-1i*aa)
    ## J <- exponential %*%data_proj-H*Msp[ch]

    J<-temp$J 
    S<-apply(Conj(J)*J,1,mean)
    list(f=f,S=Re(S))
}
