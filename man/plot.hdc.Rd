\name{plot.hdc}
\alias{plot.hdc}
\title{
Cattel's scree-test for 'hdc' class objects.
}
\description{
This function plots Cattel's scree-test of parameters coming from using hdda or hddc functions.
}
\usage{
\method{plot}{hdc}(x, threshold = 0.2, ...)
}

\arguments{
  \item{x}{
A 'hdc' class object obtained using hdda() or hddc() methods.
}
  \item{threshold}{
The threshold used in Cattel's scree-test.
}
\item{...}{
Arguments based from or to other methods.
}
}
\value{
The plot of the eigen values and of the sequential differences of the eigen values. The dimension to retain is the 
}
\note{
It is the same graph obtained using the option graph=TRUE in hdda() function.
}
\references{
C. Bouveyron, S. Girard and C. Schmid, High-Dimensional Data Clustering, Computational Statistics and Data Analysis, vol. 52 (1), pp. 502-519, 2007
C. Bouveyron, S. Girard and C. Schmid, High Dimensional Discriminant Analysis, Communications in Statistics : Theory and Methods, vol. 36 (14), pp. 2607-2623, 2007
}
\author{
R. Aidan, L. Berge, C. Bouveyron, S. Girard 
}

\seealso{
\code{\link{hdda}}, \code{\link{hddc}}, \code{\link{predict.hdc}}.
}
\examples{
# Example 1 :
data(wine)
a <- wine[,-1]
z <- wine[,1]
prms1 <- hdda(a,z,model='best',scaling=TRUE,d="B")
#the plot of Cattell's scree-test with a threshold of 0.2:
plot(prms1)     
#it shows the dimension obtained with a 0.3 threshold:
plot(prms1,0.3)                                    


}
\keyword{ cattell }
\keyword{ hdda }
\keyword{ hddc }
\keyword{ clustering }
