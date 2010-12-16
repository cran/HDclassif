pck_hdda_bic <-
function(data,model,par,mod){
	k<-par$k
	d<-par$d
	b<-par$b
	a<-par$a
	mu<-par$mu
	N<-nrow(data)
	p<-ncol(data)
	prop<-par$prop
	
	if(length(b)==1) b<-rep(b,length=k)
	if (length(a)==1) a<-matrix(a,k,max(d))
	else if (length(a)==k) a<-matrix(a,k,max(d))
	else if (model=='AJBQD') a<-matrix(a,k,d[1],byrow=TRUE)
	
	b[b<1e-10]<-1e-10
	
	if (mod){
		som_a=c()
		for (i in 1:k) som_a[i]<-sum(log(a[i,1:d[i]]))
		L<- -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop) + p*(1+log(2*pi))))
	}
	else if (any(model==c("ABQD","AJBQD"))){
		Q<-rep(list(par$Q),k)
		K<-matrix(0,k,N)
		for (i in 1:k) {
			s<-sum(log(a[i,1:d[i]]))
			X<-data-matrix(mu[i,],N,p,byrow=TRUE)
			proj<-(X%*%Q[[i]])%*%t(Q[[i]])
			A<-(-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B<-X-proj
			K[i,]<-rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
		A<--1/2*t(K)
		L<-sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))/N
	}
	else L<-par$loglik[length(par$loglik)]
	
	
	ro<-k*p+k-1
	tot<-sum(d*(p-(d+1)/2))
	D<-sum(d)
	d<-d[1]
	to<-d*(p-(d+1)/2)
	if (model=='AKJBKQKDK') m<-ro+tot+2*k+D
	else if (model=='AKBKQKDK') m<-ro+tot+3*k
	else if (model=='ABKQKDK') m<-ro+tot+2*k+1
	else if (model=='AKJBQKDK') m<-ro+tot+k+D+1
	else if (model=='AKBQKDK') m<-ro+tot+2*k+1
	else if (model=='ABQKDK') m<-ro+tot+k+2
	else if (model=='AKJBKQKD') m<-ro+k*(to+d+1)+1
	else if (model=='AKBKQKD') m<-ro+k*(to+2)+1
	else if (model=='ABKQKD') m<-ro+k*(to+1)+2
	else if (model=='AKJBQKD') m<-ro+k*(to+d)+2
	else if (model=='AKBQKD') m<-ro+k*(to+1)+2
	else if (model=='ABQKD') m<-ro+k*to+3
	else if (model=='AJBQD') m<-ro+to+d+2
	else if (model=='ABQD') m<-ro+to+3
	bic<--2*L+m*log(N)/N
	return(-bic)
}

