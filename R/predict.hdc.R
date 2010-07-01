predict.hdc <-
function(object,data,cls=NULL,...){
	p<-ncol(data)
	N<-nrow(data)
	k<-object$k
	a<-object$a
	b<-object$b
	mu<-object$mu
	d<-object$d
	prop<-object$prop
	Q<-object$Q
	x<-as.matrix(data)
	if (length(N)==0) {
		N<-1
		p<-length(data)
		x<-matrix(data,N,p)
	}
	if (length(object$scaling)!=0){
		x<-scale(x,center=object$scaling$mu,scale=object$scaling$sd)
	}
	
	if(length(b)==1) b<-rep(b,length=k)
	if (length(a)==1) a<-matrix(a,k,max(d))
	else if (length(a)==k) a<-matrix(a,k,max(d))
	else if (object$model=='AIBQD') a<-matrix(a,k,d[1],byrow=TRUE)
	
	b[b<1e-10]<-1e-10
	
	if(object$model=="AIBQD") {
		K<-diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if (object$model=="ABQD") {
		K<-diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
	}
	else{
		K<-matrix(0,k,N)
		for (i in 1:k) {
			s<-sum(log(a[i,1:d[i]]))
			X<-x-matrix(mu[i,],N,p,byrow=TRUE)
			proj<-(X%*%Q[[i]])%*%t(Q[[i]])
			A<-(-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B<-X-proj
			K[i,]<-rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])
		}
	}
	
	t<-matrix(0,N,k,dimnames=list(1:N,1:k))
	for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
	result=max.col(t)
	
	if (!is.null(object$kname)){
		result<-factor(result,labels=object$kname,levels=seq_along(object$kname))
		colnames(t)<-object$kname
	}
	if (!is.null(cls)){
		if (length(unique(result))!=length(unique(cls))) cat("Warning : there is one or more classes not represented\n")
		else pck_hdda_tclass(cls,result)
	}
	class(t)<-'hd'
	list(class=result,posterior=t)
}

