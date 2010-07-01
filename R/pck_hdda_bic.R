pck_hdda_bic <-
function(data,model,par){
	k<-par$k
	d<-par$d
	b<-par$b
	a<-par$a
	N<-nrow(data)
	p<-ncol(data)
	prop<-par$prop
	
	if(length(b)==1) b<-rep(b,length=k)
	if (length(a)==1) a<-matrix(a,k,max(d))
	else if (length(a)==k) a<-matrix(a,k,max(d))
	else if (model=='AIBQD') a<-matrix(a,k,d[1],byrow=TRUE)
	
	b[b<1e-10]<-1e-10
	
	L<-0
	for (i in 1:k){
		som_a<-sum(log(a[i,1:d[i]]))
		L<-L+prop[i]*(som_a+(p-d[i])*log(b[i])-2*log(prop[i])+p*(1+log(2*pi)))
	}	
	
	ro<-k*p+k-1
	tot<-sum(d*(p-(d+1)/2))
	D<-sum(d)
	d<-d[1]
	to<-d*(p-(d+1)/2)
	if (model=='AKIBKQKDK') m<-ro+tot+2*k+D
	else if (model=='AKBKQKDK') m<-ro+tot+3*k
	else if (model=='ABKQKDK') m<-ro+tot+2*k+1
	else if (model=='AKIBQKDK') m<-ro+tot+k+D+1
	else if (model=='AKBQKDK') m<-ro+tot+2*k+1
	else if (model=='ABQKDK') m<-ro+tot+k+2
	else if (model=='AKIBKQKD') m<-ro+k*(to+d+1)+1
	else if (model=='AKBKQKD') m<-ro+k*(to+2)+1
	else if (model=='ABKQKD') m<-ro+k*(to+1)+2
	else if (model=='AKIBQKD') m<-ro+k*(to+d)+2
	else if (model=='AKBQKD') m<-ro+k*(to+1)+2
	else if (model=='ABQKD') m<-ro+k*to+3
	else if (model=='AIBQD') m<-ro+to+d+2
	else if (model=='ABQD') m<-ro+to+3
	bic<-L+m*log(N)/N
	return(-bic)
}

