\name{predict.hdc}
\alias{predict.hdc}
\title{
Prediction method for 'hdc' class objects.
}
\description{
This function computes the class prediction of a dataset with respect to the model-based supervised and unsupervised classification methods \code{\link{hdda}} and \code{\link{hddc}}.
}
\usage{
\method{predict}{hdc}(object, data, cls = NULL, ...)
}

\arguments{
  \item{object}{
An 'hdc' class object obtained by using hdda() or hddc() function.
}
  \item{data}{
A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{cls}{
It is optional and used to be compared to the predicted classes, default is FALSE.
}
\item{...}{
Arguments based from or to other methods.
}
}
\value{
 \item{class}{vector of the predicted class.}
 \item{prob}{The matrix of the probabilities to belong to a class for each observation and each class.}
If the true class vector is given in parameters (i.e. cls!=F), it also gives the good classification rate and the confusion matrix.
}
\references{
C. Bouveyron, S. Girard and C. Schmid, High-Dimensional Data Clustering, Computational Statistics and Data Analysis, vol. 52 (1), pp. 502-519, 2007
}
\author{
R. Aidan, L. Berge, C. Bouveyron, S. Girard 
}

\seealso{
\code{\link{hdda}}, \code{\link{hddc}}.
}
\examples{
# Example 1 :
datas<-simuldata(1000,1000,50)
X<-datas$X
clx<-datas$clx
Y<-datas$Y
cly<-datas$cly
prms1=hddc(X,k=3,graph=TRUE,algo="CEM",init='param')                 #clustering of datahdlearn
prms1$class                                                          #class vector obtained by the clustering
res1=predict(prms1,X,clx)                                            #only to see the good classification rate and the confusion matrix
res2=predict(prms1,Y)         
res2$class                                                           #the class predicted using hddc parameters on the test dataset

# Example 2 :
data(Crabs)
prms3=hddc(Crabs[,-1],k=4,graph=TRUE,algo="EM",init='kmean')         #clustering of the Crabs dataset
res3=predict(prms3,Crabs[,-1],Crabs[,1])

}
\keyword{ hddc }
\keyword{ hdda }
\keyword{ clustering }