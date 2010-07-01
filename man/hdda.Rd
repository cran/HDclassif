\name{HDDA}
\alias{hdda}
\title{
High Dimensional Discriminant Analysis
}
\description{
HDDA is a model-based discriminant analysis method assuming each class of the dataset live in a proper Gaussian subspace which is much smaller than the original one, the hdda.learn function calculates the parameters of each subspace in order to predict  the class of new observation of this kind.
}
\usage{
hdda(data, cls, model = "AkiBkQkDk", graph = FALSE, d = "Cattell", 
threshold = 0.2, show = TRUE, scaling = FALSE)
}

\arguments{
  \item{data}{
A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{cls}{
The vector of the class of each observations, its type can be numeric or string.
}
  \item{model}{
The models used with the HDDA method are the following : "AkiBkQkDk", "AkBkQkDk", "ABkQkDk", "AkiBQkDk", "AkBQkDk", "ABQkDk", "AkiBkQkD", "AkBkQkD", "ABkQkD", "AkiBQkD", "AkBQkD", "ABQkD", "AiBQD", "ABQD", "best" and "dbest".

Aki are the parameters of the classes subspaces. If Aki : each class have its parameters and there is one parameter for each dimension ; if Ak : the classes have different parameters but there is only one per class ; if Ai : all the classes have the same parameters for each dimension (it's a particular case with a common orientation matrix) ; if A : all classes have the same one parameter.

Bk are the noises of the classes subspaces. If Bk : each class have its proper noise ; if B : all classes have the same noise.

Qk is the orientation matrix of each class. If Qk : all classes have its proper orientation matrix ; if Q : all classes have orientation.

Dk is the intrinsic dimension of each class. If Dk : the dimensions are free and found by Cattel's Scree-Test ; if D : the dimension is common to all classes and is either found using Cattel's Scree-Test on covariance matrix of the whole dataset, or specified by the user.

The model "best" will compute the six first models (with free dimension, Dk), give their BIC and keep the model with the highest BIC value.

The model 'dbest' will compute the eigth last models (with common dimension, D), give their BIC and keep the model with the highest BIC value.
}
  \item{graph}{
Logical, TRUE or FALSE. If TRUE, it will display the graph with the eigen values and the Cattel's Scree-Test for each class. Cattel's Scree-Test is to find the dimension of each class subspace. It may display the BIC dimension choice if it has been chosen for dimension selection (d="Bic").
}
  \item{d}{
 This parameter have several functions.
It specifies if the intrinsic dimensions of the classes are found with Cattell's scree-test or using the BIC criterion (i.e. the dimension is here found maximising the BIC criterion : the dimension kept is the one before the first fall of the BIC w.r.t. the dimension). Use d="C" for Cattell's scree-test and d="B" for using the BIC criterion.
It gives the dimension used with the models with common dimensions :"AkiBkQkD" to "ABQD", if used in this order, it must be an integer stricly inferior to the number of variables. 
Note that if a d is not specified for a common dimension model, the dimension will be gathered via the common covariance matrix, using Cattell's scree-test or the BIC criterion, wether d is equal to "Cattell" or "Bic", can be abreviated to "C" or "B".

}
  \item{threshold}{
A float stricly within 0 and 1. It is the threshold used in the Cattel's Scree-Test.
}
  \item{show}{
Use show=FALSE to settle off the informations that may be printed.
}
  \item{scaling}{
Logical : wether to scale the dataset or not.
}
}

\value{
hdda.learn returns an 'hdc' object ; it's 	a list containing :
\item{ model }{The name of the model.}
\item{ k }{The number of classes.}
\item{ d }{The dimensions of each class.}
\item{ a }{The parameters of each class subspace.}
\item{ b }{The noise of each class subspace.}
\item{ mu }{The mean of each variable for each class.}
\item{ prop }{The proportion of each class.}
\item{ ev }{The eigen values of the var/covar matrix.}
\item{ Q }{The orthogonal matrix of orientation of each class.}
\item{ kname }{The name of each class.}
\item{ BIC }{The BIC value of the model used.}
\item{ scaling }{The centers and the standard deviation of the original dataset.}
}
\references{
 C. Bouveyron, S. Girard and C. Schmid, High Dimensional Discriminant Analysis, Communications in Statistics : Theory and Methods, vol. 36 (14), pp. 2607-2623, 2007
}
\author{
R. Aidan, L. Berge, C. Bouveyron, S. Girard  
}
\seealso{
\code{\link{hddc}}, \code{\link{predict.hdc}}.
}
\examples{
#example 1 :
datas<-simuldata(1000,1000,50,k=5)
X<-datas$X
clx<-datas$clx
Y<-datas$Y
cly<-datas$cly
prms1 <- hdda(X,clx,graph=TRUE)         #we get the HDDA parameters

cl1 <- predict(prms1,Y,cly)
cl1$class                               #the class vector of Y estimated with HDDA

#example 2 :
data(wine)
a <- wine[,-1]
z <- wine[,1]
prms2 <- hdda(a,z,model='best',scaling=TRUE,d="B")
cl2 <- predict(prms2,a,z)



}
\keyword{ hdda }
\keyword{ predict }
\keyword{ classification }
