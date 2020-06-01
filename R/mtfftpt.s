mtfftpt <- function(data,tapers,nfft,t,f,findx){
#'  Multi-taper fourier transform for point process given as times
#'
#' Multi-taper fourier transform for point process given as times
#'
#' @param data
#' @param tapers from dpss
#' @param nfft  length of padded data
#' @param t  time points at which tapers are calculated
#' @param  f frequencies at which fft is caculated
#' @param  findx  index corresponding to frequencies f
#' @return J fft in form frequency index x taper index
#' @return Msp (number of spikes per sample in each channel
#' @return Nsp (number of spikes in each channel
#' 
#' @export

    
    C <-1
    K <-dim(tapers)[2]
    nfreq <- length(f); #% number of frequencies
    
    #H<-apply(tapers$v,2,fft)
    temp<-matrix(0, 1269,99)
    H<-apply(rbind(tapers,temp),2,fft)
    
    
    # dim(A) 10 20  -- apply(A,1,sum) means sum the rows. apply(A,2,sum) means sum the columns
    # However in matlab  fft(A,n,1) means apply fft along the first dimension i.e. to the columns. 
    
    
    #H <- fft(tapers$v) #% fft of tapers  fft(X,n,1); length of X is less than n, then X is padded with trailing zeros to length n.
#    dim(H) #[1] ] 2827   99
    H <- H[findx,] #; % restrict fft of tapers to required frequencies
    # 2049   99
    
    #H=fft(tapers,nfft,1);  % fft of tapers
    w <- 2*pi*f    #; % angular frequencies at which ft is to be evaluated
    sp <- 0
    Msp<- 0
    Nsp<- 0
    ch <- 1
    dtmp <- data
    
    indx <- which(dtmp>=min(t) & dtmp<=max(t))
    dtmp <- dtmp[indx]
    Nsp[ch] <- length(dtmp)
    Msp[ch]  <- Nsp[ch]/length(t);
    
    data_proj <-apply(tapers,2,fm<-function(x){pracma::interp1(t,x,dtmp) })
    
    aa<-outer(w,(dtmp-t[1]))
    exponential <- exp(-1i*aa)
    J <- exponential %*%data_proj-H*Msp[ch]
    aa<-list()
    aa$J <- J
    aa$Nsp <- Nsp
    aa$Msp <- Msp
    aa
}
