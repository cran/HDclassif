\name{HDDC}
\alias{hddc}
\title{
High Dimensional Data Clustering
}
\description{
HDDC is a model-based clustering method. It is based on the Gaussian Mixture Model and on the idea that the data lives in subspaces with a lower dimension than the dimension of the original space. It uses the Expectation - Maximisation algorithm to estimate the parameters of the model.
}
\usage{
hddc(data, K=1:10, model=c("AkjBkQkDk"), threshold=0.2, 
	com_dim=NULL, itermax=60, eps=1e-3, graph=FALSE, 
	algo='EM', d="Cattell", init='kmeans', show=TRUE,
	mini.nb=c(5,10), scaling=FALSE, min.individuals=2, noise.ctrl=1e-8, \dots)
}

\arguments{
  \item{data}{
A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{K}{
A vector of integers specifying the number of clusters for which the BIC and the parameters are to be calculated; the function keeps the parameters which maximises the BIC. Note that the length of the vector K can't be larger than 20. Default is 1:10.
}
  \item{model}{
A character string vector, or an integer vector indicating the models to be used. The available models are: "AkjBkQkDk" (default), "AkBkQkDk", "ABkQkDk", "AkjBQkDk", "AkBQkDk", "ABQkDk", "AkjBkQkD", "AkBkQkD", "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD", "AjBQD", "ABQD". It is not case sensitive and integers can be used instead of names, see details for more information. Several models can be used, if it is, only the results of the one which maximizes the BIC criterion is kept. To run all models, use model="ALL".
}
  \item{threshold}{
A float stricly within 0 and 1. It is the threshold used in the Cattell's Scree-Test.
}
	\item{com_dim}{
	It is used only for common dimensions models. The user can give the common dimension he wants. If used, it must be an integer. Its default is set to NULL.
}
  \item{itermax}{
The maximum number of iterations allowed. The default is 60.
}
  \item{eps}{
A positive double. It is the stopping criterion: the algorithm stops when the difference between two successive Log Likelihoods is lower than \sQuote{eps}.
}
  \item{graph}{
	It is for comparison sake only, when several estimations are run at the same time (either when using several models, or when using cross-validation to select the best dimension/threshold). If graph = TRUE, the plot of the results of all estimations is displayed. Default is FALSE.
}
  \item{algo}{
A character string indicating the algorithm to be used. The available algorithms are the Expectation-Maximisation ("EM"), the Classification E-M ("CEM") and the Stochastic E-M ("SEM"). The default algorithm is the "EM".
}
   \item{d}{
   Either \dQuote{Cattell} (default) or \dQuote{BIC}. See details for more information.
}
  \item{init}{
A character string or a vector of clusters. It is the way to initialize the E-M algorithm. There are five ways of initialization: \dQuote{kmeans} (default), \dQuote{param}, \dQuote{random} or \dQuote{mini-em}. See details for more information. It can also be directly initialized with a vector containing the prior classes of the observations.
}
  \item{show}{
Use show = FALSE to settle off the informations that may be printed.
}
  \item{mini.nb}{
A vector of integers of length two. This parameter is used in the \dQuote{mini-em} initialization. The first integer sets how many times the algorithm is repeated; the second sets the maximum number of iterations the algorithm will do each time. For example, if init=\dQuote{mini-em} and mini.nb=c(5,10), the algorithm wil be lauched 5 times, doing each time 10 iterations; finally the algorithm will begin with the initialization that maximizes the log-likelihood. 
}
  \item{scaling}{
Logical: whether to scale the dataset (mean=0 and standard-error=1 for each variable) or not. By default the data is not scaled.
}
  \item{min.individuals}{
This parameter is used to control for the minimum population of a class. If the population of a class becomes stricly inferior to 'min.individuals' then the algorithm stops and gives the message: 'pop<min.indiv.'. Here the meaning of "population of a class" is the sum of its posterior probabilities. The value of 'min.individuals' cannot be lower than 2.
}
	\item{noise.ctrl}{
	This parameter avoids to have a too low value of the 'noise' parameter b. It garantees that the dimension selection process do not select too many dimensions (which leads to a potential too low value of the noise parameter b). When selecting the intrinsic dimensions using Cattell's scree-test or BIC, the function doesn't use the eigenvalues inferior to noise.ctrl, so that the intrinsic dimensions selected can't be higher or equal to the order of these eigenvalues.
}
  \item{\dots}{Any arguments that can be used by the function "kmeans" when this initialization is selected (namely: algorithm; nstart; iter.max).}
}

\value{
hddc returns an 'hdc' object; it's a list containing:
\item{ model }{The name of the model.}
\item{ K }{The number of classes.}
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
\item{ likelihood }{The Log Likelihood for each iteration.}
\item{ class }{The class vector obtained by the clustering.}
\item{ posterior }{The matrix of the probabilities to belong to a class for each observation and each class.}
}

\details{
Some information on the signification of the model names:
\describe{
	\item{Akj are the parameters of the classes subspaces:}{
		\itemize{
		\item{if Akj: each class has its parameters and there is one parameter for each dimension}
		\item{if Ak: the classes have different parameters but there is only one per class}
		\item{if Aj: all the classes have the same parameters for each dimension (it's a particular case with a common orientation matrix)}
		\item{if A: all classes have the same one parameter}
		}
	}

	\item{Bk are the noises of the classes subspaces:}{
		\itemize{
			\item{If Bk: each class has its proper noise}
			\item{if B:  all classes have the same noise}
		}
	}

	\item{Qk is the orientation matrix of each class:}{ 
		\itemize{
			\item{if Qk: all classes have its proper orientation matrix}
			\item{if Q: all classes have the same orientation matrix}
		}
	}
	
	\item{Dk is the intrinsic dimension of each class:}{ 
		\itemize{
			\item{if Dk: the dimensions are free and proper to each class}
			\item{if D: the dimension is common to all classes}
		}
	}
}
The model \dQuote{ALL} will compute all the models, give their BIC and keep the model with the highest BIC value.
Instead of writing the model names, they can also be specified using an integer.  1 represents the most general model (\dQuote{AkjBkQkDk}) while 14 is the most constrained (\dQuote{ABQD}), the others  number/name matching are given below. Note also that several models can be run at once, by using a vector of models (e.g. model = c("AKBKQKD","AKJBQKDK","AJBQD") is equivalent to model = c(8,4,13); to run the 6 first models, use model=1:6). If all the models are to be run, model="all" is faster than model=1:14. 
\tabular{lcclc}{
AkjBkQkDk \tab   1   \tab   \tab  AkjBkQkD \tab   7   \cr 
AkBkQkDk \tab   2   \tab \tab  AkBkQkD \tab   8   \cr   
ABkQkDk \tab   3   \tab  \tab ABkQkD \tab   9   \cr   
AkjBQkDk \tab   4   \tab  \tab  AkjBQkD \tab   10   \cr   
AkBQkDk \tab   5   \tab  \tab  AkBQkD \tab   11   \cr   
ABQkDk \tab   6   \tab  \tab  ABQkD \tab   12  \cr
AjBQD \tab 13 \tab  \tab ABQD \tab 14
}

The parameter d, is used to select the intrinsic dimensions of the subclasses. Here are his definictions:
		\itemize{
			\item{\dQuote{Cattell}:}{
				The Cattell's scree-test is used to gather the intrinsic dimension of each class. If the model is of common dimension (models 7 to 14), the scree-test is done on the covariance matrix of the whole dataset.
			}
			\item{\dQuote{BIC}:}{
				The intrinsic dimensions are selected with the BIC criterion. See Bouveyron \emph{et al.} (2010) for a discussion of this topic.
				For common dimension models, the procedure is done on the covariance matrix of the whole dataset.
			}
			\item{Note that "Cattell" (resp. "BIC") can be abreviated to "C" (resp. "B") and that this argument is not case sensitive.}
		}
		
The different initializations are:
\describe{
	\item{\dQuote{param}:}{it is initialized with the parameters, the means being generated by a multivariate normal distribution and the covariance matrix being common to the whole sample}

	\item{\dQuote{mini-em}:}{it is an initialization strategy, the classes are randomly initialized and the E-M algorithm makes several iterations, this action is repetead a few times (the default is 5 iterations and 10 times), at the end, the initialization choosen is the one which maximise the log-likelihood (see mini.nb for more information about its parametrization)}

	\item{\dQuote{random}:}{the classes are randomly given using a multinomial distribution}

	\item{\dQuote{kmeans}:}{the classes are initialized using the kmeans function (with: algorithm="Hartigan-Wong"; nstart=4; iter.max=50); note that the user can use his own arguments for kmeans using the dot-dot-dot argument } 
	
	\item{A prior class vector:}{It can also be directly initialized with a vector containing the prior classes of the observations}
}

}



\references{
Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519

Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, \bold{46}(6), 1--29, url: \href{http://www.jstatsoft.org/v46/i06/}{http://www.jstatsoft.org/v46/i06/}
}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard
}

\seealso{
\code{\link{hdda}}, \code{\link{predict.hdc}}, \code{\link{plot.hdc}}.
}
\examples{
# Example 1:
data <- simuldata(1000, 1000, 50)
X <- data$X
clx <- data$clx
Y <- data$Y
cly <- data$cly

#clustering of the simulated dataset:
prms1 <- hddc(X, K=3, algo="CEM", init='param')                

#class vector obtained by the clustering:
prms1$class                

#We can look at the adjusted rand index to assess the goodness of fit
res1 <- predict(prms1, X, clx)                                            
res2 <- predict(prms1, Y)       
#the class predicted using hddc parameters on the test dataset:  
res2$class                                                           


# Example 2:
data(Crabs)

#clustering of the Crabs dataset:
prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='mini-em')        
res3 <- predict(prms3, Crabs[,-1], Crabs[,1])

#another example using the Crabs dataset
prms4 <- hddc(Crabs[,-1], K=1:8, graph=TRUE, model=c(1,2,7,9))

#model=c(1,2,7,9) is equivalent to model=c("AKJBKQKDK","AKBKQKDK","AKJBKQKD","ABKQKD") 
res4 <- predict(prms4, Crabs[,-1], Crabs[,1])
}
\keyword{ hddc }
\keyword{ clustering }
