pck_hddc_e_step <-
function(x,par){
	p<-ncol(x)
	N<-nrow(x)
	k<-par$k
	a<-par$a
	b<-par$b
	mu<-par$mu
	d<-par$d
	prop<-par$prop
	Q<-par$Q
	
	b[b<1e-10]<-1e-10

	if(par$model=="AIBQD") {
		K<-diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if(par$model=="ABQD") {
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
	
	t<-matrix(0,N,k)
	for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
	t
}

