#'
#' run the specificed function and save result to a file 
#' 
#' @param Fun the function name will be carried out
#' @param file a requied parameter that specify the name of the file and its path
#' @param redo a logical flag to rerun the given function no matter whether result file is existed or not
#'
#'@details
#' all of the values returned by the given function was saved under the \code{result} variable in the file with the given name.
#'
#' @return no return value
#'
#'@examples
#'
#' runAndSave(length,file="./test_runAndSave.RData",redo=TRUE,x=list(1:4,1:3))
#' 
#' file.exists("./test_runAndSave.RData") # expect true
#' 
#' load("./test_runAndSave.RData")
#' 
#' result==2 #expect true 
#' 
#' #remove the file after try the examples
#' file.remove("./test_runAndSave.RData")
#' 

runAndSave<-function(Fun,file=NULL,redo=FALSE,...){
  if(is.null(file))
    stop("a file name including path should be given")
  if(redo | !file.exists(file)){
    result=Fun(...)
    save(result,file=file)
  }
}