% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwsem.R
\name{get.AIC}
\alias{get.AIC}
\title{Title get.AIC}
\usage{
get.AIC(sem.model, MAG, data)
}
\arguments{
\item{sem.model}{A list containing the structural equations, each created
using either the gam or the gamm functions of the mgcv package}

\item{MAG}{A matrix encoding the directed acyclic graph (DAG), or the mixed
acyclic graph, of the structural equations model.  This is created using
the DAG or the makeGM function in the gmm library}

\item{data}{A data frame holding the observed data used in the calls to
the models in the sem.model object}
}
\value{
A data frame containing the log-likelihood of the full SEM (LL),
the number of free parameters that were estimated (K), along with the
AIC and the bias-corrected AIC (AICc)
}
\description{
Title get.AIC
}
\examples{
library(mgcv)
library(ggm)
set.seed(10)
N<-1000
L1<-rnorm(N)
x1<-0.5*L1+rnorm(N,0,sqrt(1-0.5^2))
x2<-0.5*x1+rnorm(N,0,sqrt(1-0.5^2))
x3<-0.5*L1+0.5*x2+rnorm(N,0,sqrt(1-2*0.5^2))
x4<-0.5*x3+rnorm(N,0,sqrt(1-0.5^2))
my.dat<-data.frame(x1,x2,x3,x4)

my.list<-list(gam(x1~1,data=my.dat),
             gam(x2~x1,data=my.dat),
             gam(x3~x2,data=my.dat),
             gam(x4~x3,data=my.dat))

full.dag<-DAG(x1~L1,x2~x1,x3~x2+L1,x4~x3)
my.mag<-DAG.to.MAG.in.pwSEM(full.DAG=full.dag,latents=c("L1"))
get.AIC(sem.model=my.list,MAG=my.mag,data=my.dat)
}
