Eloglik <- function(theta, delta1, delta2, E.gamma, E.log.gamma, w = NULL)
{
  if (!is.null(w))
  {
  out <- (1/theta)*(log(1/theta))  + (1/theta - 1) * mean(w*E.log.gamma) -
    (1 / theta) * mean(w*E.gamma) - log(gamma(1/theta))
  } else {
    out <- (1/theta)*(log(1/theta))  + (1/theta - 1) * mean(E.log.gamma) -
      (1 / theta) * mean(E.gamma) - log(gamma(1/theta))
  }
  #cat("theta = ", theta, ",   g(theta) = ", out, "\n")
  return(out)
}

# g_fun <- function(the,mean.Eg,mean.Elg) {
#   temp <- 1/the*log(1/the) + (1/the - 1)*mean.Elg -1/the*mean.Eg - log(gamma(1/the))
#   return(-temp)
# }
