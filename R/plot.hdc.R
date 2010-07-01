plot.hdc <-
function(x,threshold=0.2,...){
	ev<-x$ev
	k<-x$k
	if (is.matrix(ev)==FALSE) ev<-matrix(ev,k,length(ev),byrow=TRUE)
	X<-abs(t(diff(t(ev))))
	nbis<-ncol(X)
	d<-maxi<-Nmax<-c()
	for (i in 1:k) Nmax[i]<-max(sum(ev[i,]>0)-1,2)
	for (i in 1:k) for (j in Nmax[i]:1) {
		maxi[i]<-max(X[i,1:Nmax[i]])
		if (X[i,j]>=threshold*maxi[i]) {
			d[i]<-j
			break
		} 
		else if(j==1) d[i]<-1
	}

	par(mfrow=c(k*(k<=4)+4*(k>4),2*(1*(k%%4!=0)+floor(k/4))))
	for (i in 1:k){
		sub1<-paste("Class #",i,", d",i,"=",d[i],sep="")
		plot(ev[i,1:(min(d[i]+10,nbis))],type="h",col="green",main=paste("Ordered Eigen Values\nClass #",i,sep=""),xlab="",ylab="",lwd=3)
		plot(X[i,1:(min(d[i]+10,Nmax[i]))],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab="",xlab=paste("threshold=",threshold,sep=''))
		abline(h=threshold*maxi[i],lty=3)	
		points(d[i],X[i,d[i]],col='red')
	}
}

