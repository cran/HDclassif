simuldata <-
function(nlearn,ntest,p,K=3,prop=NULL,d=NULL,a=NULL,b=NULL){
	N=nlearn+ntest
	if (length(prop)==0) prop<-rep(1/K,K)
	else if (length(prop)!=K) stop("Proportions don't fit with the number of classes.")
	else prop<-prop/sum(prop)
	
	# Class sizes
	n<-floor(prop*N)
	N<-sum(n)

	#MEANS
	mu<-matrix(0,K,p)
	j<-sample(p,K)
	mu[cbind(1:K,j)]<-10
	
	# Intrinsic dimensions
	if ( length(d)==0 )	d<-sort(ceiling(runif(K,0,12*(p>20)+5*(p<=20 && p>=6)+(p<6)*(p-1))),decreasing=TRUE)
	else if ( length(d)!=K || !any(is.numeric(d)) ) stop("Wrong value of d.")
	
	# Orientation matrices
	Q<-vector(mode='list',length=K)
	for (i in 1:K) Q[[i]]<-qr.Q(qr(mvrnorm(p,mu=rep(0,p),Sigma=diag(1,p))))
	
	# Variance in the class-specific subspace
	if ( length(a)==0 ) a<-sort(ceiling(runif(K,30,350)))
	else if ( length(a)!=K || !any(is.numeric(a)) ) stop("Wrong value of a.")
	if ( length(b)==0 )b<-sort(ceiling(runif(K,0,25)))
	else if ( length(b)!=K || !any(is.numeric(b)) ) stop("Wrong value of b.")
	
	# Simulation
	S<-vector(mode='list',length=K)
	for (i in 1:K)	S[[i]]<-crossprod(Q[[i]]%*%sqrt(diag(c(rep(a[i],d[i]),rep(b[i],p-d[i])))))
	
	cls<-X<-NULL
	for (i in 1:K)	X<-rbind(X,mvrnorm(n[i],mu=mu[i,],Sigma=S[[i]]))
	
	for (i in 1:K) cls<-c(cls,rep(i,n[i]))
	
	ind<-sample(1:N,N)
	prms<-list(a=a,b=b,prop=prop,d=d,mu=mu)
	data <- list(X=X[ind[1:nlearn],],clx=cls[ind[1:nlearn]],Y=X[ind[(nlearn+1):N],],cly=cls[ind[(nlearn+1):N]],prms=prms)
	
}

