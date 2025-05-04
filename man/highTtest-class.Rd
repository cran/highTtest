\name{highTtest-class}
\Rdversion{1.1}
\docType{class}
\alias{highTtest-class}
\alias{BH,highTtest-method}
\alias{CK,highTtest-method}
\alias{pi_alt,highTtest-method}
\alias{pvalue,highTtest-method}
\alias{ST,highTtest-method}

\title{Class \code{"highTtest"}}
\description{
Value object returned by call to \code{highTtest()}. 
}
\section{Objects from the Class}{
This object should not be created by users.
}
\section{Slots}{
  \describe{
    \item{\code{CK}:}{Object of class \code{matrix} or NULL.
                      A matrix of logical values. The
                      rows correspond to features, ordered as
                      provided in input \code{dataSet1}. The columns correspond to
                      levels of significance. Matrix elements are TRUE if
                      feature was determined to be significant 
                      by the Cao-Kosorok method.
                      The significance value associated with each column is
                      dictated by the input \code{gammas}.}
    \item{\code{pi1}:}{Object of class \code{numeric} or NULL. 
                       The estimated proportion of alternative hypotheses
                       calculated using the Cao-Kosorok method.}
    \item{\code{pvalue}:}{Object of class \code{numeric}. 
                       The vector of p-values calculated using the
                       two-sample t-statistic. }
    \item{\code{ST}:}{Object of class \code{matrix} or NULL. 
                      If requested, a matrix of logical values. The
                      rows correspond to features, ordered as
                      provided in input \code{dataSet1}. The columns correspond to
                      levels of significance. Matrix elements are TRUE if
                      feature was determined to be significant 
                      by the Storey-Tibshirani (2003) method.
                      The significance value associated with each column is
                      dictated by the input \code{gammas}.}
    \item{\code{BH}:}{Object of class \code{matrix} or NULL 
                      If requested, A matrix of logical values. The
                      rows correspond to features, ordered as
                      provided in input \code{dataSet1}. The columns correspond to
                      levels of significance. Matrix elements are TRUE if
                      feature was determined to be significant 
                      by the Benjamini-Hochberg (1995) method.
                      The significance value associated with each column is
                      dictated by the input \code{gammas}.}
    \item{\code{gammas}:}{Object of class \code{numeric}.  
                      Vector of significant values provided as 
                      input for the calculation.}
  }
}
\section{Methods}{
  \describe{
    \item{BH}{\code{signature(x = "highTtest")}: 
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Benjamini-Hochberg
              (1995) method. }
    \item{CK}{\code{signature(x = "highTtest")}: 
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Cao-Kosorok
              (2011) method. }
    \item{pi_alt}{\code{signature(x = "highTtest")}: Retrieves the
        estimated proportion of alternative hypotheses
        obtained by the Cao-Kosorok (2011) method. }
    \item{plot}{\code{signature(x = "highTtest")}: Generates a plot
        of the number of significant features as a function of the
        level of significance as calculated for each method (CK,BH, and/or
         ST) }
    \item{pvalue}{\code{signature(x = "highTtest")}: Retrieves the
        vector of p-values calculated using the
                       two-sample t-statistic. }
    \item{ST}{\code{signature(x = "highTtest")}:  
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Storey-Tibshirani
             (2003) method. }
    \item{vennD}{\code{signature(x = "highTtest")}: Generates 
        two- and three-dimensional Venn diagrams comparing the
        features selected by each method. Implements methods of
        package \pkg{VennDiagram}. In addition to the \code{highTtest}
        object, the level of significance, \code{gamma}, must
        also be provided.
	 }
}}
\references{
Cao, H. and Kosorok, M. R. (2011).
Simultaneous critical values for t-tests in very high dimensions.
{Bernoulli}, 17, 347--394. PMCID: PMC3092179.

Benjamini, Y. and Hochberg, Y. (1995).
Controlling the false discovery rate: A practical and powerful approach to 
multiple testing.
{Journal of the Royal Statistical Society: Series B}, 57, 289--300.

Storey, J. and Tibshirani, R. (2003).
Statistical significance for genomewide studies.
{Proceedings of the National Academy of Sciences, USA}, 100, 9440--9445.


}
\author{
Authors: Hongyuan Cao, Michael R. Kosorok, and 
Shannon T. Holloway <shannon.t.holloway@gmail.com> 
Maintainer: Shannon T. Holloway <shannon.t.holloway@gmail.com> 
}

\examples{
showClass("highTtest")
}
\keyword{classes}
