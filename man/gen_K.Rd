\name{gen_K}
\alias{gen_K}
\title{gen_K}

\description{Generate a sample of a Gaussian Mixture Model whose centers are generate randomly on a sphere of radius \code{radius}.
}
\usage{
gen_K(n=500,d=5,K=3,pcont=0,df=1,
      cont="Student",min=-5,max=5,radius=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{A positive integer giving the number of data per cluster. Default is \code{500}.}
\item{d}{A positive integer giving the dimension. Default is \code{5}.}
\item{K}{A positive integer giving the number of clusters. Default is \code{3}.}
\item{pcont}{A scalar between \code{0} and \code{1} giving the proportion of contaminated data.}
\item{df}{A positive integer giving the degrees of freedom of the law of the contaminated data if \code{cont='Student'}. Default is \code{1}.}
\item{cont}{The law of the contaminated data. Can be \code{'Student'} (default) and \code{'Unif'}.}
\item{min}{A scalar giving the lower bound of the uniform law if \code{cont='Unif'}. Default is \code{-5}. }
\item{max}{A scalar giving the upper bound of the uniform law if \code{cont='Unif'}. Default is \code{5}.}
\item{radius}{The radius of the sphere on each the centers of the class are generated. Default is \code{5}.}
}
\value{
A list with:
\item{X}{A numerical matrix giving the generated data.}
\item{cluster}{An character vector specifying the true classification.}
}

\seealso{
See also \code{\link{Kmedians}} and \code{\link{Kmeans}}.
}
\keyword{Gaussian Mixture Model}
\examples{
n <- 500
K <- 3
pcont <- 0.2
ech <- gen_K(n=n,K=K,pcont=pcont)
X=ech$X
}
