#'
#'The relative root mean square error
#'
#'@param the real value
#'@param estimates a vector of estimates of the real value
#'

RRMSE=function(real,estimates){
  return(sqrt(mean((real-estimates)^2))/real)
}
