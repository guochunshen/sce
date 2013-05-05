\name{knq}
\alias{knq}
\title{To find the kth nearest neighborhood in the q equal angle sector for each forcal point/event}
\usage{
  knq(focals, obs, k = 1, q = 1, type = "ptoe")
}
\arguments{
  \item{focals}{a data.frame with colnames x and y that
  contains spaital locaiton of each focals}

  \item{obs}{the observed population distribution, a
  data.frame with colnames x and y}

  \item{k}{the order of neighborhood}

  \item{q}{the number of equal angle sectors}

  \item{type}{the type of neighborhood distance. it can be
  "ptoe" (point to event distance) or "etoe" (event to
  event distance)}
}
\description{
  To find the kth nearest neighborhood in the q equal angle
  sector for each forcal point/event
}

