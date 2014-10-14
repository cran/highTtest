\name{highTtest}
\alias{highTtest}
\alias{plot,highTtest}
\alias{vennD,highTtest}
\alias{CK,highTtest}
\alias{ST,highTtest}
\alias{BH,highTtest}
\alias{pi_alt,highTtest}
\alias{pvalue,highTtest}
\alias{CK}
\alias{ST}
\alias{BH}
\alias{pi_alt}
\alias{pvalue}

\title{
Simultaneous critical values for t-tests in very high dimensions}

\description{
Implements the method developed by Cao and Kosorok (2011) for the significance 
analysis of thousands of features in high-dimensional biological studies. 
It is an asymptotically valid data-driven procedure to find critical values 
for rejection regions controlling the k-familywise error rate, false discovery 
rate, and the tail probability of false discovery proportion.}

\usage{
highTtest(dataSet1, dataSet2, gammas, compare = "BOTH", cSequence = NULL, 
tSequence = NULL)
}

\arguments{
  \item{dataSet1 }{data.frame or matrix containing the dataset for subset 1
                   for the two-sample t-test.}
  \item{dataSet2 }{data.frame or matrix containing the dataset for subset 2
                   for the two-sample t-test.}
  \item{gammas   }{vector of significance levels at which feature 
                   significance is to be determined.}
  \item{compare  }{one of ("ST", "BH", "Both", "None"). 
                   In addition to the Cao-Kosorok method, 
                   obtain feature significance indicators using
                   the Storey-Tibshirani method (ST) 
                   (Storey and Tibshirani, 2003),
                   the Benjamini-Hochberg method (BH), 
                   (Benjamini andHochberg, 1995),
                   "both" the ST and the BH methods, or 
                   do not consider alternative methods (none).}
  \item{cSequence}{A vector specifying the values of c to be considered
                   in estimating the proportion of alternative hypotheses. 
                   If no vector is provided, a default of 
                   seq(0.01,6,0.01) is used. See Section 2.3 of
                   Cao and Kosorok (2011) for more information.}
  \item{tSequence}{A vector specifying the search space for the critical 
                   t value. If no vector is provided, a default of 
                   seq(0.01,6,0.01) is used.}
}
\details{
The Storey-Tibshirani (2003), ST, method implemented in \pkg{highTtest} is
adapted from the implementation written by Alan Dabney and John D. Storey and 
available from 

http://www.bioconductor.org/packages/release/bioc/html/qvalue.html. 

The comparison capability is included only for convenience and reproducibility 
of the original manuscript. For a complete analysis based on the ST method, 
the user is referred to the \pkg{qvalue} package available through the
bioconductor archive.

The following methods retrieve individual results from a highTtest
object, x: 

    \code{BH(x)}: 
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Benjamini-Hochberg
              (1995) method.

    \code{CK(x)}: 
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Cao-Kosorok
              (2011) method.

    \code{pi_alt(x)}: Retrieves the
        estimated proportion of alternative hypotheses
        obtained by the Cao-Kosorok (2011) method. 

    \code{pvalue(x)}: Retrieves the
        vector of p-values calculated using the
                       two-sample t-statistic. 

    \code{ST(x)}:  
              Retrieves a matrix of logical values. The
              rows correspond to features, the columns to levels
              of significance. Matrix elements are TRUE if feature
              was determined to be significant by the Storey-Tibshirani
             (2003) method. 



A simple x-y plot comparing the number of significant
features as a function of the level significance level
can be generated using

    \code{plot(x,...)}: Generates a plot
        of the number of significant features as a function of the
        level of significance as calculated for each method (CK,BH, and/or
         ST). Additional plot controls can be passed through the ellipsis.

When comparisons to the ST and BH methods are requested, Venn diagrams
can be generated using 

    \code{vennD(x, gamma, ...)}: Generates 
        two- and three-dimensional Venn diagrams comparing the
        features selected by each method. Implements methods of
        package \pkg{colorfulVennPlot}. In addition to the \code{highTtest}
        object, the level of significance, \code{gamma}, must
        also be provided. Most control argument of the
        \pkg{colorfulVennPlot} package can be passed through the ellipsis.


}
\value{
Returns an object of class \code{highTtest}.

}
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
Shannon T. Holloway <sthollow@ncsu.edu> 
Maintainer: Shannon T. Holloway <sthollow@ncsu.edu> 
}

\examples{
set.seed(123)
x1 <- matrix(c(runif(500),runif(500,0.25,1)),nrow=100)
obj <- highTtest(dataSet1=x1[,1:5], 
                 dataSet2=x1[,6:10], 
                 gammas=seq(0.1,1,0.1),
                 tSequence=seq(0.001,3,0.001))

#Print number of significant features identified in each method
colSums(CK(obj))
colSums(ST(obj))
colSums(BH(obj))

#Plot the number of significant features identified in each method
plot(obj, main="Example plot")

vennD(obj, 0.8, Title="Example vennD")

#Proportion of alternative hypotheses
pi_alt(obj)

#p-values
pvalue(obj)

}

