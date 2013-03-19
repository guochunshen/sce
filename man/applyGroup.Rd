\name{applyGroup}
\alias{applyGroup}
\title{apply the same function on every member of a group, multicore is supported on Linux with multicores CPU}
\usage{
  applyGroup(com, group, Fun, verbpro = FALSE,
    mc.cores = 5, multicore = FALSE, repository = NULL,
    ...)
}
\arguments{
  \item{com}{a \code{\link{scp}} object.}

  \item{group}{factor, defination of unit of individuals
  used to fit model. e.g. species name, DBH class}

  \item{Fun}{a function name that will apply on every
  element of a group}

  \item{verbpro}{logical flage to show the verbose of
  process of calculation}

  \item{mc.cores}{number of cores used in parallel
  calculation on Linux}

  \item{...}{other parameters passed to the \code{Fun}
  function}

  \item{repository}{an repository object that can save
  middle result into the repository file.  see
  \code{\link{create_repository}} for detail. it only
  usefull to run long term calculation if the calculation
  environmental is not stable and will crash sometime.}
}
\value{
  a list of results that returned by the function
}
\description{
  apply the same function on every member of a group,
  multicore is supported on Linux with multicores CPU
}

