\name{multiplot}
\alias{multiplot}
\title{Multiple plot function}
\usage{
  multiplot(..., plotlist = NULL, file, cols = 1,
    layout = NULL)
}
\arguments{
  \item{...}{ggplot objects can be passed in ..., or to
  plotlist (as a list of ggplot objects)}

  \item{cols}{Number of columns in layout}

  \item{layout}{A matrix specifying the layout. If present,
  'cols' is ignored.}
}
\description{
  Multiple plot function
}
\details{
  If the layout is something like matrix(c(1,2,3,3),
  nrow=2, byrow=TRUE), then plot 1 will go in the upper
  left, 2 will go in the upper right, and 3 will go all the
  way across the bottom.
}

