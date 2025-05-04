\name{vennD-methods}
\docType{methods}
\alias{vennD}
\alias{vennD-methods}
\alias{vennD,highTtest,numeric-method}
\title{ ~~ Methods for Function \code{vennD}  ~~}
\description{
Generates 2- or 3-dimensional Venn diagrams comparing the
features selected by the Cao-Kosorok method to those selected
by the Storey-Tibshirani (2003) method 
and/or the Benjamini-Hoshberg (1995) method.
This S4 method is simply a wrapper 
for \code{draw.pairwise.venn()} and \code{draw.triple.venn()} of 
package \pkg{VennDiagram}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "highTtest", gamma="numeric", ...)}}{
Object returned by a call to \code{highTtest()}.
\code{gamma} is the level of significance. Additional control variables
for the methods of \code{draw.pairwise.venn()} and \code{draw.triple.venn()} of 
package \pkg{VennDiagram} can be passed through the ellipsis.
}
}}
\keyword{methods}

