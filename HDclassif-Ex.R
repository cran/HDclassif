pkgname <- "HDclassif"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('HDclassif')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("hdda")
### * hdda

flush(stderr()); flush(stdout())

### Name: HDDA
### Title: High Dimensional Discriminant Analysis
### Aliases: hdda
### Keywords: hdda predict classification

### ** Examples

#example 1:
data<-simuldata(1000, 1000, 50, K=5)
X <- data$X
clx <- data$clx
Y <- data$Y
cly <- data$cly
#we get the HDDA parameters:
prms1 <- hdda(X, clx)         

cl1 <- predict(prms1, Y, cly)
#the class vector of Y estimated with HDDA:
cl1$class                     

#another model is used:
prms1 <- hdda(X, clx, model=12)
#model=12 is equivalent to model="ABQkD"     
cl1 <- predict(prms1, Y, cly) 

#example 2:
data(wine)
a <- wine[,-1]
z <- wine[,1]
prms2 <- hdda(a, z, model='all', scaling=TRUE, d="B", graph=TRUE)
cl2 <- predict(prms2, a, z)

#getting the best dimension
#using a common dimension model
#we do LOO-CV using cv.vfold=nrow(a)
prms3 <- hdda(a, z, model="akjbkqkd", d="CV", cv.vfold=nrow(a), scaling=TRUE, graph=TRUE)

cl3 <- predict(prms3, a, z)




cleanEx()
nameEx("hddc")
### * hddc

flush(stderr()); flush(stdout())

### Name: HDDC
### Title: High Dimensional Data Clustering
### Aliases: hddc
### Keywords: hddc clustering

### ** Examples

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

#only to see the good classification rate and the confusion matrix:                                          
res1 <- predict(prms1, X, clx)                                            
res2 <- predict(prms1, Y)       
#the class predicted using hddc parameters on the test dataset:  
res2$class                                                           


# Example 2:
data(Crabs)

#clustering of the Crabs dataset:
prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='param')        
res3 <- predict(prms3, Crabs[,-1], Crabs[,1])

#another example on the Crabs dataset
prms4 <- hddc(Crabs[,-1], K=1:8, graph=TRUE, model=c(1,2,7,9))

#model=c(1,2,7,9) is equivalent to model=c("AKJBKQKDK","AKBKQKDK","AKJBKQKD","ABKQKD") 
res4 <- predict(prms4, Crabs[,-1], Crabs[,1])



cleanEx()
nameEx("plot.hdc")
### * plot.hdc

flush(stderr()); flush(stdout())

### Name: plot.hdc
### Title: Cattell's Scree-Test for 'hdc' class objects.
### Aliases: plot.hdc
### Keywords: cattell hdda hddc clustering

### ** Examples

# Example 1 :
data(wine)
a <- wine[,-1]
z <- wine[,1]

prms1 <- hdda(a, z, model="AkBkQkDk", scaling=TRUE, d="B")

#the plot related to the selection that has been done: BIC
plot(prms1)     

#it shows the plot of Cattell's scree-test, with a threshold of .3
plot(prms1,"Cattell",0.3)                         


prms2 <- hdda(a, z, model="AkBkQkD", scaling=TRUE, d="c")
#the plot related to the selection that has been done: Cattell's scree-test
plot(prms2) 
#the plot of the BIC
plot(prms2,"b") 




cleanEx()
nameEx("predict.hdc")
### * predict.hdc

flush(stderr()); flush(stdout())

### Name: predict.hdc
### Title: Prediction method for 'hdc' class objects.
### Aliases: predict.hdc
### Keywords: hddc hdda clustering

### ** Examples

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

#only to see the good classification rate and the confusion matrix:                     
res1 <- predict(prms1, X, clx)                                            
res2 <- predict(prms1, Y)       

#the class predicted using hddc parameters on the test dataset:  
res2$class                                                           


# Example 2:
data(Crabs)
#clustering of the Crabs dataset:
prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='kmeans')        
res3 <- predict(prms3, Crabs[,-1], Crabs[,1])

#the confusion matrix:
res3$confusion




cleanEx()
nameEx("simuldata")
### * simuldata

flush(stderr()); flush(stdout())

### Name: simuldata
### Title: Gaussian Data Generation
### Aliases: simuldata
### Keywords: generation gaussian

### ** Examples

data <- simuldata(500, 1000, 50, K=5, prop=c(0.2,0.25,0.25,0.15,0.15))
X <- data$X
clx <- data$clx
f <- hdda(X, clx)
Y <- data$Y
cly <- data$cly
e <- predict(f, Y, cly)




cleanEx()
nameEx("wine")
### * wine

flush(stderr()); flush(stdout())

### Name: wine
### Title: Wine dataset
### Aliases: wine
### Keywords: datasets

### ** Examples

data(wine)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
