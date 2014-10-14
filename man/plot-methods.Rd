\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,highTtest-method}
\title{ ~~ Methods for Function \code{plot}  ~~}
\description{
 Generates a simple x-y plot giving the number of significant
features as a function of the level of significance. If comparisons
to Storey-Tibshirani and Benjamini-Hochberg methods were requested
by the user, these will automatically be included in the plot. 
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Plot method as implemented by other packages.
}

\item{\code{signature(x = "highTtest")}}{
Object returned by a call to \code{highTtest()}.
}
}}
\keyword{methods}

