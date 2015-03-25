\name{MPG-package}
\alias{MPG-package}
\alias{MPG}
\docType{package}
\title{
Comparison across mixture models through mixtures of perturbed Gaussians
}
\description{
Bayesian model for comparison across mixtures of perturbed Gaussians.
The model detects variations in the mixture weights and in the mixture locations across the mixture models. 
}
\details{
\tabular{ll}{
Package: \tab MPG\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-10\cr
License: \tab GPL (>= 3)\cr
}

}
\author{
Jacopo Soriano 

Maintainer: Jacopo Soriano <jacopo.soriano@gmail.com>
}
\references{
Literature or other references for background information
}

\keyword{ package }
\seealso{
 Optional links to other man pages, e.g.
 \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
}
\examples{

n = c(250, 250)
p = 4

Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))

ans = mpg(Y, C)  
plotDiff(ans, type = "weight")
plotDiff(ans, type = "shift")
}
