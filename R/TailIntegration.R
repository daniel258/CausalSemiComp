#################################################################
## Function for CausalSemiComp
### S integration of a function
# This is generalization of KM integration for any tail function S
##################################################################

## S integration
# f- function,
# times - time points in which S was estimated
# S - tail function
# tau - optional positive value for RMST-like calcualtions
# start - time zero
# f.start - relevant only if f is non-null, the value of f at time zero
TailIntegration <- function(times = NULL, S, S.points = T, tau = NULL, f = NULL, f.points = T, start = 0, f.start = 0)
{
  if(is.null(tau)) {
    warning(paste0("End-time taken to be the last time point, setting tau = ", times[length(times)]),
            call. = F)
    tau <- times[length(times)]
  }
  if(!(length(unique(times))==length(times))) stop("times should not include duplicated items")
  if(is.unsorted(times,strictly = T)) stop("times should be sorted in an increasing order")
   if(tau > times[length(times)]) {
     time.pts <- c(times, tau)
     if (S.points==T) {
       S.at.pts <- c(S, 0) # the zero is not used for the calculation
     } else {
       S.at.pts <- c(sapply(times, S), 0)
     }
     if (!is.null(f)) {
      if (f.points==T) {
         f.at.pts <- c(f, 0) # the zero is not used for the calculation
      } else {
         f.at.pts <- c(sapply(times, f), 0)
      }}
   } else if (tau==times[length(times)]) {
     time.pts <- times
     if (S.points==T) {
       S.at.pts <- S
     } else {
       S.at.pts <- sapply(time.pts, S)
     }
     if (!is.null(f)) {
       if (f.points==T) {
        f.at.pts <- f
        } else {
        f.at.pts <- sapply(time.pts, f)
     }}
   } else {
     times <- times[times <= tau]
     if(times[length(times)] != tau)
     {
       time.pts <- c(times, tau)
     } else {
       time.pts <- times
     }
     if (S.points==T) {
       S.at.pts <- S[1:length(time.pts)]
     } else {
       S.at.pts <- sapply(time.pts, S)
     }
     if (!is.null(f)) {
       if (f.points==T) {
         f.at.pts <- f[1:length(time.pts)]
       } else {
         f.at.pts <- sapply(time.pts, f)
      }}
   }
  time.diff <- diff(c(start, time.pts))
  if (is.null(f))
  {
    rectang <- time.diff * c(1, S.at.pts[-length(S.at.pts)])
  } else {
    rectang <- time.diff * c(f.start, f.at.pts[-length(S.at.pts)]) * c(1, S.at.pts[-length(S.at.pts)])
  }
  res <- sum(rectang)
  return(res)
}
