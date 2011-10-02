#Clustering example on the "Crabs" dataset :
#Here two models are fitted and different number of clusters 
#are tested. The model and the number of clusters that 
#maximize the BIC criterion are selected.

data(Crabs)
X<-Crabs[,-1]
clx<-Crabs[,1]
prms<-hddc(X,k=1:6,model=c("AkjBkQkDk","ABkQkD"))
prms
res<-predict(prms,X,clx)
prms<-hddc(X,k=2:6,model=c(4:5,10),cgraph=TRUE)

#Now a PCA on the "Crabs" dataset is done on the two first 
#principal axis.
#It is an illustration of the clustering process. 
#The means of the clusters are represented by the points 
#and the directions of each group is represented by a line. 
#All the steps of the algorithm are shown.
#The algorithm, the initialization and the model can be chosen.
demo_hddc()




