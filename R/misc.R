#' @description
#' Functions to retrieve angles and coordinates
#' 
get_coord <- function(theta){
  return(c(cos(theta), sin(theta)))
}

get_theta <- function(omega){
  return( atan2(omega[2], omega[1]) )
}