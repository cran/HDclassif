plot.hdc <-
function(x,threshold=0.2,...){
	ev<-x$ev
	k<-x$k
	p<-ncol(ev)
	if (is.matrix(ev)==FALSE) ev<-matrix(ev,k,length(ev),byrow=TRUE)
	dvp<-abs(t(diff(t(ev))))
	nbis<-ncol(dvp)
	d<-maxi<-Nmax<-c()
	for (i in 1:k) Nmax[i]<-max(sum(ev[i,]>0)-1,2)
	if (p==2) d=rep(1,k)
	else{
		for (i in 1:k){
			Nmax[i]<-max(p-2,2)
			maxi[i]<-threshold*max(dvp[i,1:Nmax[i]])
			d[i]<-max(which(dvp[i,1:Nmax[i]]>=maxi[i]))
		}
	}
	par(mfrow=c(k*(k<=4)+4*(k>4)-(any(k==c(5,6,10))),2*(1*(k%%4!=0)+floor(k/4))))
	for (i in 1:k){
		sub1<-paste("Class #",i,", d",i,"=",d[i],sep="")
		plot(ev[i,1:(min(d[i]+10,nbis))],type="h",col="green",main=paste("Ordered Eigen Values\nClass #",i,sep=""),xlab="",ylab="",lwd=3)
		plot(dvp[i,1:(min(d[i]+10,Nmax[i]))],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab="",xlab=paste("threshold=",threshold,sep=''))
		abline(h=maxi[i],lty=3)	
		points(d[i],dvp[i,d[i]],col='red')
	}
}

