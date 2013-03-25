#'
#'divide a large plot into small quadrats
#'
#'@param com a scp object
#'@param nx number of quadrat in the x axis
#'@param ny number of qudart in the y axis
#'
#'@return
#'a new scp object with a "ploti" trait in the traits of the community
#'it also has an attribute "spaced" contains spatial distance between quadrats
#'
#'@examples
#'data(testData)
#'
#'newdata=quadratize(testData,5,5)
#'
#' table(newdata$traits$ploti)
#'


quadratize=function(com,nx,ny){
  data.ppp=com$com
  win=data.ppp$win
  xi_breaks=seq(0-1e-5,win$xrange[2]+1e-5,length.out=nx+1)
  yi_breaks=seq(0-1e-5,win$yrange[2]+1e-5,length.out=ny+1)
  xi=as.numeric(cut(data.ppp$x,breaks=xi_breaks))
  yi=as.numeric(cut(data.ppp$y,breaks=yi_breaks))
  ploti=xi+(yi-1)*nx
  com$traits$ploti=ploti
  #spatial distance between quadrat
  plotx=rep((xi_breaks[-1]+xi_breaks[-(nx+1)])/2,times=ny)
  ploty=rep((yi_breaks[-1]+yi_breaks[-(nx+1)])/2,each=nx)
  plotd=dist(cbind(plotx,ploty))
  attr(com,"spaced")=plotd
  return(com)
}
