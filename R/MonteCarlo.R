MonteCarlo <-
function(p,M=10^4)
{
  if(M%%1 !=0 | M<1) stop("M must be a positive integer")
  if(length(p)<1) stop("p length must greater or equal than 1")
  if(max(p) > 1 | min(p) < 0) stop("p vector must contain real numbers in [0,1]")

  cat("   ","\n")
  cat("###########################################","\n")
  cat(" Opportunistic Model - Monte Carlo results ","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Number of Hops N:  ",length(p)),"\n")
  cat(paste("Probabilities:     ", paste(p, collapse=", ")),"\n")
  cat(paste("Simulations M:     ", paste(M, collapse=", ")),"\n")
  cat("   ","\n")
  out = MC(p,M)
  cat("   ","\n")
  print(out)
  return(out)
}
