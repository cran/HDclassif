\name{predict.hdc}
\alias{predict.hdc}
\title{
Prediction method for \sQuote{hdc} class objects.
}
\description{
This function computes the class prediction of a dataset with respect to the model-based supervised and unsupervised classification methods \code{\link{hdda}} and \code{\link{hddc}}.
}
\usage{
\method{predict}{hdc}(object, data, cls = NULL, \dots)
}

\arguments{
  \item{object}{
An 'hdc' class object obtained by using hdda or hddc function.
}
  \item{data}{
A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{cls}{
It is optional and used to be compared to the predicted classes, default is NULL.
}
\item{\dots}{
Arguments based from or to other methods. Not currently used.
}
}
\value{
	\item{class}{vector of the predicted class.}
	\item{prob}{The matrix of the probabilities to belong to a class for each observation and each class.}
	If the initial class vector is given to the argument \sQuote{cls} then the adjusted rand index (ARI) is also returned. Also the following object is returned:
	\item{ARI}{The confusion matrix of the classification.}
}
\references{
Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607--2623

Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519

Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, \bold{46}(6), 1--29, url: \href{http://www.jstatsoft.org/v46/i06/}{http://www.jstatsoft.org/v46/i06/}
}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard
}

\seealso{
\code{\link{hdda}}, \code{\link{hddc}}.
}
\examples{
# Example 1:
data <- simuldata(1000, 1000, 50)
X <- data$X
clx <- data$clx
Y <- data$Y
cly <- data$cly

#clustering of the gaussian dataset:
prms1 <- hddc(X, K=3, algo="CEM", init='param')      
           
#class vector obtained by the clustering:
prms1$class                   

# only to see the good classification rate and 
# the Adjusted Rand Index:                     
res1 <- predict(prms1, X, clx)                                            
res2 <- predict(prms1, Y)       

#the class predicted using hddc parameters on the test dataset:  
res2$class                                                           


# Example 2:
data(Crabs)
#clustering of the Crabs dataset:
prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='kmeans')        
res3 <- predict(prms3, Crabs[,-1], Crabs[,1])

}
\keyword{ hddc }
\keyword{ hdda }
\keyword{ clustering }
