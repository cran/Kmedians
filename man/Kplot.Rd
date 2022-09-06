\name{Kplot}
\alias{Kplot}
\title{Kplot}

\description{A plot function for K-medians and K-means
}
\usage{
Kplot(a,propplot=0.95,graph=c('Two_Dim','Capushe','Profiles','SE','Criterion'),
      bestresult=TRUE,Ksel=FALSE,bycluster=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{a}{Output from \code{\link{Kmedians}} or \code{\link{Kmeans}}.}
\item{propplot}{A scalar between \code{0} and \code{1} giving the propotion of data considered for the different graphs.}
\item{graph}{A string specifying the type of graph requested.
Default is \code{c('Two_Dim','Capushe','Profiles','SE','Criterion')}.}
\item{bestresult}{A logical indicating if the graphs must be done for the result chosen by the selected criterion. Default is \code{TRUE}.}
\item{Ksel}{A logical or positive integer giving the chosen number of clusters for each the graphs should be drawn.}
\item{bycluster}{A logical indicating if the data selected for \code{'Two_Dim'} and \code{'Profiles'} graphs should be selected by cluster or not. Default is \code{TRUE}.}
}

\value{No return value.}

\seealso{
See also \code{\link{Kmedians}} and \code{\link{Kmeans}}.
}
\keyword{Robust clustering}
\examples{
n <- 500
K <- 3
pcont <- 0.2
ech <- gen_K(n=n,K=K,pcont=pcont)
X <-ech$X
res <- Kmedians(X,par=FALSE)
Kplot(res)
}
