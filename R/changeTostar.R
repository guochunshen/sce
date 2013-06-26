#' A help functions to change pvalue into the correspondence significant level using *
#' 
#' @param x the pvalue
#'
#' @examples
#'changeTostar(1)
#'
#'changeTostar(0.06)
#'
#'changeTostar(0.03)
#'
#'changeTostar(0.001)


changeTostar<-function(x){
  threestars=x<0.01
  twostars=x>=0.01 & x<0.05
  onestar=x>=0.05 & x<0.1
  zerostar=x>=0.1
  x[threestars]="***"
  x[twostars]="**"
  x[onestar]="*"
  x[zerostar]=""
  return(x)
}