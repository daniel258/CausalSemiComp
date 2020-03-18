# Auxiliry function for Beran curves. Replace NA in a vector x with minimum value
RepNAmin <- function(x)
{
  min.x <- min(x, na.rm = T)
  x[is.na(x)] <- min.x
  return(x)
}
