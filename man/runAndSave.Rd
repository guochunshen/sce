\name{runAndSave}
\alias{runAndSave}
\title{run the specificed function and save result to a file}
\usage{
  runAndSave(Fun, file = NULL, redo = FALSE, ...)
}
\arguments{
  \item{Fun}{the function name will be carried out}

  \item{file}{a requied parameter that specify the name of
  the file and its path}

  \item{redo}{a logical flag to rerun the given function no
  matter whether result file is existed or not}
}
\value{
  no return value
}
\description{
  run the specificed function and save result to a file
}
\details{
  all of the values returned by the given function was
  saved under the \code{result} variable in the file with
  the given name.
}
\examples{
runAndSave(length,file="./test_runAndSave.RData",redo=TRUE,x=list(1:4,1:3))

file.exists("./test_runAndSave.RData") # expect true

load("./test_runAndSave.RData")

result==2 #expect true

#remove the file after try the examples
file.remove("./test_runAndSave.RData")
}

