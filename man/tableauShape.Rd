% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tableaux.R
\name{tableauShape}
\alias{tableauShape}
\title{Shape of a tableau}
\usage{
tableauShape(tableau)
}
\arguments{
\item{tableau}{a tableau (list of vectors having the same mode)}
}
\value{
The shape of the tableau. This is an integer partition whose 
  \code{i}-th part is the number of boxes in the \code{i}-th row of 
  the tableau.
}
\description{
The shape of a tableau.
}
\examples{
tableau <- list(c(2, 1, 3), c(5, 2))
tableauShape(tableau)
}
