mtspectrumpt<-function(PP, Fs  =   3000,err  =  c(2, 0.0500),fpass  =  c(0, 1500),
                       pad  =  0,tt =  matrix(c(50, 99),1,2),tapers=matrix(c(50, 99),1,2),fscorr =  0){
    mintime <- min(PP)
    maxtime <- max(PP)
    dt <- 1/Fs #% sampling time
    t <- seq(from=mintime-dt,by=dt,to=maxtime+dt)
    
    N<-length(t)# % number of points in grid for dpss
    nfft <- max(2^(pracma::nextpow2(N)+pad),N) #number of points in fft of prolates
    
    aa<-getfgrid(Fs,nfft,fpass)
    f<-aa$f
    findx<-aa$findx
    
    tapers <- dpsschk(tapers,N,Fs) #check tapers 
#    sz<- dim( tt)
#    if (sz[1]==1 && sz[2]==2){
#        tapers <-multitaper::dpss(N,tt[2],tt[1])
#        tapers$v = tapers$v*sqrt(Fs)
#    }
    
    #[J,Msp,Nsp]=mtfftpt(data,tapers,nfft,t,f,findx); % mt fft for point process times
    C <-1
    K <-dim(tapers)[2]
    nfreq <- length(f); #% number of frequencies
    
    #H<-apply(tapers$v,2,fft)
    temp<-matrix(0, 1269,99)
    H<-apply(rbind(tapers,temp),2,fft)
    
    
    # dim(A) 10 20  -- apply(A,1,sum) means sum the rows. apply(A,2,sum) means sum the columns
    # However in matlab  fft(A,n,1) means apply fft along the first dimension i.e. to the columns. 
    
    
    #H <- fft(tapers$v) #% fft of tapers  fft(X,n,1); length of X is less than n, then X is padded with trailing zeros to length n.
    dim(H) #[1] ] 2827   99
    H <- H[findx,] #; % restrict fft of tapers to required frequencies
    # 2049   99
    
    #H=fft(tapers,nfft,1);  % fft of tapers
    w <- 2*pi*f    #; % angular frequencies at which ft is to be evaluated
    sp <- 0
    Msp<- 0
    Nsp<- 0
    ch <- 1
    dtmp <- S1
    
    indx <- which(dtmp>=min(t) & dtmp<=max(t))
    dtmp <- dtmp[indx]
    Nsp[ch] <- length(dtmp)
    Msp[ch]  <- Nsp[ch]/length(t);
    
    data_proj <-apply(tapers,2,fm<-function(x){pracma::interp1(t,x,dtmp) })
    
    aa<-outer(w,(dtmp-t[1]))
    exponential <- exp(-1i*aa)
    J <- exponential %*%data_proj-H*Msp[ch]
    S<-apply(Conj(J)*J,1,mean)
}
