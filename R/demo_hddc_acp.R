demo_hddc_acp <-
function(X,z,hd=NULL,...){
	Q=hd$Q
	MU=hd$mu
	d=hd$d
	if (is.matrix(Q)) {
		svg=vector(mode='list',length=4)
		for (i in 1:4) svg[[i]]=as.matrix(Q)
		Q=svg
	}

	X=as.matrix(X)
	p=ncol(X)
	k=max(z)
	
	CO=cov(X)
	Z=-eigen(CO,symmetric=T)$vectors
	coul=1:4
	
	V=X%*%Z
	patch=c(3,4,8,20)
	plot(V[,1],V[,2],col=z,pch=patch[z],...)
	

	proj=matrix(,k,p)
	for (i in 1:k) proj[i,]=tcrossprod(Q[[i]])%*%matrix(10,p,1)+MU[i,]
	
	x=proj%*%Z
	y=MU%*%Z
	points(y[,1],y[,2],col=coul[1:k],pch=19,lwd=7)
	
	for (i in 1:k) {
		pente=(x[i,2]-y[i,2])/(x[i,1]-y[i,1])
		oo=x[i,2]-pente*x[i,1]
		xb=(2*y[i,1]-sqrt(50^2/(pente^2+1)))/2
		xa=(2*y[i,1]+sqrt(50^2/(pente^2+1)))/2
		lines(c(xa,xb),oo+pente*c(xa,xb),col=coul[i],type='l')
	}

}

