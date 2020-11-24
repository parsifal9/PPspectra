#' simulate events from a homogeneous Poisson process
#'
#' 
#'
#' @param rate
#' @param num.events
#' 
#' @return aa1 vector of event times
#' 
#' @export
hpp.event.times<-function (rate, num.events, t0 = 0,tT=1,sig=4)
{
    x = matrix(rexp(n = num.events , rate = rate),  num.events)
    aa1 <- t0 + apply(x, 2, cumsum)
    aa1<-round(aa1,sig)
    aa1<-aa1[aa1< tT]
    aa1<-unique(aa1)
    aa1
}
