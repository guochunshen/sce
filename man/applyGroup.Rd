\name{applyGroup}
\alias{applyGroup}
\title{apply the same function on every member of a group}
\usage{
  applyGroup(com, group, Fun, verbpro = TRUE, ...)
}
\arguments{
  \item{com}{a \code{\link{scp}} object.}

  \item{group}{factor, defination of unit of individuals
  used to fit model. e.g. species name, DBH class}

  \item{Fun}{a function name that will apply on every
  element of a group}

  \item{verbpro}{logical flage to show the verbose of
  process of calculation}

  \item{...}{other parameters passed to the \code{Fun}
  function}
}
\value{
  a list of results that returned by the function
}
\description{
  apply the same function on every member of a group
}

