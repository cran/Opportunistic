\name{routes}
\alias{routes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Routing distribution for an Opportunistic network
}
\description{
It provides the different
possible routes, their frequency as well as their respective
probabilities when considering uncertain probabilities lying
on a interval \code{p} +- \code{delta}.
}
\usage{
routes(p, delta = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
vector of probabilities of length \code{N} where \code{N} represents the number of hops
}
  \item{delta}{
Delta value when considering uncertain probabilities. The interval is of the type \code{p} \eqn{+-} \code{delta}.
}
}
\details{
By default, \code{delta} is considered to be zero disregarding uncertainty.
}
\value{
A data frame containing the routes, frequencies, and respective probabilities.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{Expected}}, \code{\link{MonteCarlo}}
}
\examples{
\dontrun{
#An N=7 Opportunistic system with probabilities p1 = 0.7,...,p7 = 0.1

> p = seq(0.7,0.1,length.out = 7)
> routes(p)
            Freq Probability   Value
route 1        1        p1^7 0.08235
route 2        6     p1^5*p2 0.10084
route 3       10   p1^3*p2^2 0.12348
route 4        4     p1*p2^3  0.1512
route 5        5     p1^4*p3 0.12005
route 6       12  p1^2*p2*p3   0.147
route 7        3     p2^2*p3    0.18
route 8        3     p1*p3^2   0.175
route 9        4     p1^3*p4  0.1372
route 10       6    p1*p2*p4   0.168
route 11       2       p3*p4     0.2
route 12       3     p1^2*p5   0.147
route 13       2       p2*p5    0.18
route 14       2       p1*p6    0.14
route 15       1          p7     0.1
Total         64
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Opportunistic}
\keyword{transmissions}
\keyword{routing}
\keyword{receptions}
\keyword{boradcast}
\keyword{network}
