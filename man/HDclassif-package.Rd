\name{HDclassif-package}
\alias{HDclassif-package}
\alias{HDclassif}
\docType{package}
\title{
High Dimensional Discriminant Analysis and Data Clustering
}
\description{
New disciminant analysis and data clustering methods for high dimensional data, based on the asumption that high-dimensional data live in different subspaces with low dimensionality, proposing a new parametrization of the Gaussian mixture model which combines the ideas of dimension reduction and constraints on the model.
}
\details{
\tabular{ll}{
Package: \tab HDclassif\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2011-07-30\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package is used to make efficient supervised and unsupervised classification with high dimensional data.
The supervised method uses the \var{hdda} function to get the data parameters and the \var{predict} function to realise the class prediction of a dataset.
The unsupervised method is implemented in the \var{hddc} function, and once the parameters are estimated, the \var{predict} gives the class prediction of other datasets. The method used in the hddc is based on the Expectation - Maximisation algorithm.
}
\author{
L.Berge C.Bouveyron S.Girard 

Maintainer: L.Berge <laurent.berge at hotmail.fr> 
}
\references{
Bouveyron C, Girard S and Schmid C, dQuote{High Dimensional Discriminant Analysis, Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607-2623, 2007

Bouveyron C, Girard S and Schmid C, \dQuote{High-Dimensional Data Clustering, Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502-519, 2007
}

\keyword{ package }