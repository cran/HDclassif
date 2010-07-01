pck_hddc_dim <-
function(DATA,k,threshold){
	p<-ncol(DATA)
	N<-nrow(DATA)
	MU<-colMeans(DATA)
	
	if (N<p) {
		Y<-(DATA-matrix(MU,N,p,byrow=TRUE))/sqrt(N)
		YYt<-tcrossprod(Y)
		ev<-eigen(YYt,symmetric=TRUE,only.values=TRUE)$values
	}
	else{
		S<-crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
		ev<-eigen(S,symmetric=TRUE,only.values=TRUE)$values
	}
	
	x<-abs(diff(ev))
	Nmax<-max(min(floor(N/k-2),p-2),2)
	for (j in Nmax:1) {
		if (x[j]>=threshold*max(x[1:Nmax])) {
			d<-j
			break
		} 
		else if(j==1) d<-1
	}
	
	if (d>min(N,p)-1) d<-min(N,p)-1
	return(d)
}

