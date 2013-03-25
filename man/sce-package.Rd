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
    V. \tab hypothesis tests and model diagnostics \cr
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
  
  \bold{To simulate a spatial population/community pattern}
  \tabular{ll}{
    \code{\link{rmypoispp}} \tab simulate a null zero point pattern given by the lambda and also avoiding too extreme point pattern \cr
    \code{\link{rCluster}} \tab simulate a cluster population distribution from a fitted cluster model or given parameters \cr
    \code{\link{rCommunity}} \tab simulate a community uniquely driven by random/dispersal/habitat/compeition \cr
  }
   \bold{Typical spatial community data within the package}
  
  \tabular{ll}{
    \code{\link{NATrees}} \tab the maps of North American trees \cr
    \code{\link{testData}} \tab a small community data for fast testing \cr
    \code{\link{BCI}} \tab the 2005 full individual mapped data from BCI plot \cr
    
  } 
  
  \bold{To manipulate a spatial population/community pattern}
  
  \tabular{ll}{
    \code{\link{subset.scp}} \tab extract some part of the community \cr
    \code{\link{removeSpecies}} \tab remove species in the given species list from the community \cr
    \code{\link{removeRareSpecies}} \tab remove rare species from the community \cr
    \code{\link{selSpeciesByAbund}} \tab select a community with specific species abundance \cr
    \code{\link{applyGroup}} \tab apply a function on each member of a group devided from a community, paralla computation supported \cr
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
    \code{\link{updateCluster}} \tab update a given cluster model manually or automatically. \cr
  }
}

\section{V. HYPOTHESIS TESTS AND MODEL DIAGNOSTICS}{

  \tabular{ll}{
    \code{\link{gofTest}} \tab goodness-of-fit test for a fitted cluster model \cr
    \code{\link{sigHabitatTest}} \tab testing significance of habitat variables by considering aggregative residual \cr
    \code{\link{sigAggreResidualTest}} \tab testing significance of aggregation in the residual. \cr
    \code{\link{backwardStep}} \tab select a best clustering model for a group of individuals based on pvalue. \cr
    \code{\link{varDecomp}} \tab calculate proportion of variance explained by habitat, internal clustering and random process. \cr
    \code{\link{envelopeVar}} \tab calculate the confidence interval of proportion of variance \cr
  }
}

\section{VI. Other functions}{

  \tabular{ll}{
    \code{\link{create_repository}} \tab save result in each step into a file for avoding unexpected error disturbance. it becomes very useful and can save you lots of time in long term calculation \cr
    \code{\link{onePlotVarDecompRoutine}} \tab runing all nessary calcution of variance decompstion for each species in one plot \cr
    \code{\link{runAndSave}} \tab run the specificed function and save result to a file \cr
    \code{\link{dataframeToIms}} \tab convert a data.frame into a list of im objects \cr
    \code{\link{downloadNATree}} \tab download the north American tree distribution files from interent \cr
    \code{\link{NAtreeRaster}} \tab transfer the zipped shapfile distribution of tree data into a dataframe of species occurance \cr
  }
}
\author{
Guochun Shen, Jian Zhang

Maintainer: Guochun Shen <shenguochun@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ spatial community ecology, spatial variance decomposition }

