pck_hdda_prms_dim <-
function(DATA,cls,threshold,graph,method){
	p<-ncol(DATA)
	N<-nrow(DATA)
	n<-prop<-c()
	k<-max(cls)
	n<-as.vector(table(cls))
	prop<-n/N

	mu<-matrix(,k,p)
	for (i in 1:k) mu[i,]<-colMeans(DATA[cls==i,])
	
	if (N<p) {
		Y<-matrix(0,N,p)
		for (i in 1:k) Y[which(cls==i),]<-(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(N)
		YYt<-tcrossprod(Y)
		ev<-eigen(YYt,symmetric=TRUE,only.values=TRUE)$values
	}
	else{
		W<-matrix(0,p,p)
		for (i in 1:k) W<-W+prop[i]*crossprod(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i]
		ev<-eigen(W,symmetric=TRUE,only.values=TRUE)$values
	}
	
	if (method=='C') {
		x<-abs(diff(ev))
		Nmax<-max(min(floor(N/k-2),p-2),2)
		for (j in Nmax:1) {
			if (x[j]>=threshold*max(x[1:Nmax])) {
				d<-j
				break
			} 
			else if(j==1) d<-1
		}
		if (graph) {
			nbis<-length(ev)-1
			x<-abs(diff(ev))
			x11()
			par(mfrow=c(2,1))
			plot(ev[1:(min(d+10,nbis))],type="h",col="green",main="Ordered Eigen Values of the \nwhole dataset covariance matrix",xlab="",ylab="",lwd=3)
			plot(x[1:(min(d+10,nbis-1))],type="l",col="blue",main=paste("Cattell's Scree-Test\nd=",d,sep=''),xlab=paste("Threshold=",threshold,sep=''),ylab='')
			abline(h=threshold*max(x),lty=3)	
			points(d,x[d],col='red')
		}
	}
	else {
		d<-0
		ev[ev<1e-10]<-1e-10

		maxi<-min(min(n[i],p))-1
		B<-c()
		for (kdim in 1:maxi){
			if (d!=0 & kdim>d+10) break
			a<-sum(ev[1:kdim])/kdim
			b<-sum(ev[(kdim+1):min(N,p)])/(p-kdim)
			if (b<1e-10) b<-1e-10
			L2<--1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))
			B[kdim]<--2*L2+(p+kdim*(p-(kdim+1)/2)+1)*log(N)/N
			if ( (d==0 & kdim>1) && B[kdim-1]<B[kdim] )	d<-kdim-1
		}
		if (d==0) d<-min(n)-1
		if (graph){
			x11()
			plot(-B,type='l',col=4,main=paste("BIC evolution w.r.t. the dimension\nd=",d,sep=''),ylab='BIC',xlab="Dimension")
			points(d,-B[d],col=2)
		}	
	}
	if (d>min(n)-1) d<-min(n)-1
	return(d)
}

