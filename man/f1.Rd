% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f1.R
\name{f1}
\alias{f1}
\title{utility function}
\usage{
f1(aa, bb, cc)
}
\arguments{
\item{aa}{vector}

\item{bb}{vector}

\item{cc}{matrix}
}
\value{
J matrix
}
\description{
returns  (exp(-1i* aa %*% t(bb)))%*% cc, The same as  exp(-1i*outer(aa,bb))%*% cc but the outer product aa %*% t(bb) may be very large
}
