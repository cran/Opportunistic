\name{Expected}
\alias{Expected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical broadcasts/transmissions/receptions for an Opportunistic model
}
\description{
This function computes the probability of success and the expected values of the number of broadcasts, transmissions and receptions for an Opportunistic model.
}
\usage{
Expected(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
vector of probabilities of length \code{N} where \code{N} represents the number of hops
}
}
\value{
A matrix with the probabilities and expected values for an Opportunistic model for all hops sizes <= N
}
\references{
Biswas, S., & Morris, R. (2004). Opportunistic routing in multi-hop wireless networks. ACM SIGCOMM Computer Communication Review, 34(1), 69-74.
}
\author{
Christian E. Galarza and Jonathan M. Olate
}
\seealso{
\code{\link{routes}}, \code{\link{MonteCarlo}}
}
\examples{
#An N=3 Opportunistic system with probabilities p = c(0.0,0.4,0.1)
res1 = Expected(p=c(0.9,0.4,0.1))
res1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Opportunistic}
\keyword{transmissions}
\keyword{receptions}
\keyword{boradcast}
\keyword{network}
