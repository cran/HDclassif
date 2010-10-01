#Supervised classification of the "wine" dataset.
#The data is scaled using the command scaling=TRUE.
#The learning is done using a random sample of 100 individuals
#whereas the testing is done with the 78 remaining individuals.
#The graph of the choice of the intrinsic dimensions is shown.
data(wine)
X<-wine[,-1]
clx<-wine[,1]
ind<-sample(178,100)
prms<-hdda(X[ind,],clx[ind],model="best",scaling=TRUE,graph=TRUE)
res<-predict(prms,X[-ind,],clx[-ind])