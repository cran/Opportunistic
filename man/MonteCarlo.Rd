\name{MonteCarlo}
\alias{MonteCarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo broadcasts/transmissions/receptions for an Opportunistic model
}
\description{
This function estimates via Monte Carlo the probability of success and the expected values of the number of broadcasts, transmissions and receptions for an Opportunistic model.
}
\usage{
MonteCarlo(p, M = 10^4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
vector of probabilities of length \code{N} where \code{N} represents the number of hops
}
  \item{M}{
Total number of Monte Carlo simulations
}
}
\details{
N is computed from \code{p} length. M is code{10^4} by default.
}
\value{
A vector with the success probability and expected values (broadcast, transmissions and receptions) for an N Opportunistic model.
}
\references{
Biswas, S., & Morris, R. (2004). Opportunistic routing in multi-hop wireless networks. ACM SIGCOMM Computer Communication Review, 34(1), 69-74.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{routes}}, \code{\link{Expected}}
}
\examples{
#Monte Carlo simulation for an N=3 Opportunistic system with probabilities
#p = c(0.0,0.4,0.1)

res2 = MonteCarlo(p=c(0.9,0.4,0.1),M=10^4)
res2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Opportunistic}
\keyword{transmissions}
\keyword{receptions}
\keyword{boradcast}
\keyword{network}
