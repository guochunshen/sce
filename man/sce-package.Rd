\name{sce-package}
\alias{sce-package}
\alias{sce}
\docType{package}
\title{
The \pkg{sce} package for spatial community ecology
}
\description{
This is a summary of the features of \pkg{sce}, a package in \code{R} for the spatial explicit anaylsis of community ecology.
}
\details{
Currently, it deals mainly with populations and community data type.

The package supports
\itemize{
    \item creation, manipulation and plotting of patterns of populations and communities
    \item exploratory data analysis 
    \item simulation of populations under various null models
    \item parametric model-fitting
    \item hypothesis tests and model diagnostics
  }
}
\section{Functions and Datasets}{
Following is a summary of the main functions and datasets in the \pkg{sce} package. Alternatively an alphabetical list of all functions and datasets is available by typing \code{library(help=sce)}.

For further information on any of these functions, type \code{help(name)} where name is the name of the function or dataset.
}

\section{Contents:}{
\tabular{ll}{
    I. \tab Creating and manipulating population/community data \cr
    II. \tab Exploratory Data Analysis \cr
    III. \tab Model fitting \cr
    IV. \tab Simulation \cr 
  }
}

\section{I. CREATING AND MANIPULATING DATA}{

  \bold{Types of data in spatial community ecology}
  
   The main types of spatial data supported by \pkg{sce} are:
   
   \tabular{ll}{
    \code{\link{scp}} \tab spatial community pattern \cr
    }
    
  \bold{To create a spatial population/community pattern}
  
  \tabular{ll}{
    \code{\link{scp}} \tab create a spatial community pattern \cr
  }
  
  \bold{To manipulate a spatial population/community pattern}
  
  \tabular{ll}{
    \code{\link{removeRareSpecies}} \tab remove rare species from the community \cr
  }
  
}

\section{II. EXPLORATORY DATA ANALYSIS }{

  \tabular{ll}{
  \code{\link{pcf_adaptive}} \tab estimate the pair correlation function by an adaptive method. \cr
  }
}


\section{III. MODEL FITTING}{
  
  \tabular{ll}{
    \code{\link{fitCluster}} \tab fitting cluster model for a spatial population/community pattern \cr
  }
}

\author{
Guochun Shen

Maintainer: Guochun Shen <shenguochun@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
}
