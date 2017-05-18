Expected <-
function(p)
{
  if(length(p)<1) stop("p length must greater or equal than 1")
  if(max(p) > 1 | min(p) < 0) stop("p vector must contain real numbers in [0,1]")

  cat("   ","\n")
  cat("###########################################","\n")
  cat(" Opportunistic Model - Theoretical results ","\n")
  cat("###########################################","\n")
  cat("   ","\n")

  cat(paste("Number of Hops N:  ",length(p)),"\n")
  cat(paste("Probabilities:     ", paste(p, collapse=", ")),"\n")
  cat("   ","\n")
  out2 = expected(p)
  cat("   ","\n")
  print(out2)
  cat("   ","\n")
  return(out2)
}
