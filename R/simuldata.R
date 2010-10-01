simuldata <-
function(nlearn,ntest,p,k=3,prop=NULL,d=NULL,a=NULL,b=NULL){
	N=nlearn+ntest
	if (length(prop)==0) prop<-rep(1/k,k)
	else if (length(prop)!=k) stop("Proportions don't fit with the number of classes.")
	else prop<-prop/sum(prop)
	
	# Class sizes
	n<-floor(prop*N)
	N<-sum(n)

	#MEANS
	mu<-matrix(0,k,p)
	j<-sample(p,k)
	mu[cbind(1:k,j)]<-10
	
	# Intrinsic dimensions
	if ( length(d)==0 )	d<-sort(ceiling(runif(k,0,12*(p>20)+5*(p<=20 && p>=6)+(p<6)*(p-1))),decreasing=TRUE)
	else if ( length(d)!=k || !any(is.numeric(d)) ) stop("Wrong value of d.")
	
	# Orientation matrices
	Q<-vector(mode='list',length=k)
	for (i in 1:k) Q[[i]]<-qr.Q(qr(mvrnorm(p,mu=rep(0,p),Sigma=diag(1,p))))
	
	# Variance in the class-specific subspace
	if ( length(a)==0 ) a<-sort(ceiling(runif(k,30,350)))
	else if ( length(a)!=k || !any(is.numeric(a)) ) stop("Wrong value of a.")
	if ( length(b)==0 )b<-sort(ceiling(runif(k,0,25)))
	else if ( length(b)!=k || !any(is.numeric(b)) ) stop("Wrong value of b.")
	
	# Simulation
	S<-vector(mode='list',length=k)
	for (i in 1:k)	S[[i]]<-crossprod(Q[[i]]%*%sqrt(diag(c(rep(a[i],d[i]),rep(b[i],p-d[i])))))
	
	cls<-X<-NULL
	for (i in 1:k)	X<-rbind(X,mvrnorm(n[i],mu=mu[i,],Sigma=S[[i]]))
	
	for (i in 1:k) cls<-c(cls,rep(i,n[i]))
	
	ind<-sample(1:N,N)
	prms<-list(a=a,b=b,prop=prop,d=d,mu=mu)
	data <- list(X=X[ind[1:nlearn],],clx=cls[ind[1:nlearn]],Y=X[ind[(nlearn+1):N],],cly=cls[ind[(nlearn+1):N]],prms=prms)
	
}

