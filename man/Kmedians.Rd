\name{Kmedians}
\alias{Kmedians}
\title{Kmedians}

\description{K-medians algorithms.
}
\usage{
Kmedians(X,nclust=1:15,ninit=0,niter=20,
         method='Offline', init=TRUE,par=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{A numerical matrix giving the data.}
\item{nclust}{A vector of positive integers giving the possible numbers of clusters. Default is \code{1:15}.}
\item{ninit}{A non negative integer giving the number of random initializations. Default is \code{0}. }
\item{niter}{A positive integer giving the number of iterations for the EM algorirthms. Default is \code{20}.}
\item{method}{The selected method for the K-medians algorithm. Can be \code{'Offline'} (default), \code{'Semi-Online'} or \code{'Online'}.}
\item{init}{A logical argument telling if the function \code{'genie'} is used for initializing the algorithm. Default is \code{TRUE}. }
\item{par}{A logical argument telling if the parallelization of the algorithm is allowed. Default is \code{TRUE}.}
}
\value{
A list with:
\item{bestresults}{A list giving all the results for the clustering selected by \code{'capushe'}.}
\item{allresults}{A list containing all the results.}
\item{SE}{A vector giving the Sum of Errors for each considered number of clusters.}
\item{cap}{The results given by the function \code{'capushe'} if \code{nclust} is  of length larger than \code{10}.}
\item{Ksel}{An integer giving the number of clusters selected by \code{'capushe'} if \code{nclust} is  of length larger than \code{10}.}
\item{data}{A numerical matrix giving the data.}
\item{nclust}{A vector of positive integers giving the considered numbers of clusters.}
For the lists \code{bestresult} and \code{allresults}:
\item{cluster}{A vector of positive integers giving the clustering.}
\item{centers}{A numerical matrix giving the centers of the clusteres.}
\item{SE}{An integer giving the Sum of Errors.}
}

\seealso{
See also \code{\link{Kmeans}}, \code{\link{Kplot}} and \code{\link{gen_K}}.
}
\references{Godichon-Baggioni, A. and Surendran, S. A penalized criterion for selecting the number of clusters for K-medians. \emph{arxiv.org/abs/2209.03597}
}
\keyword{Robust clustering}
\examples{
\dontrun{
n <- 500
K <- 3
pcont <- 0.2
ech <- gen_K(n=n,K=K,pcont=pcont)
X <-ech$X
res <- Kmedians(X,par=FALSE)
Kplot(res)
}
}
