predict.hdc  <- 
function(object,data,cls=NULL,...){
	p <- ncol(data)
	N <- nrow(data)
	K <- object$K
	a <- object$a
	b <- object$b
	mu <- object$mu
	d <- object$d
	prop <- object$prop
	Q <- object$Q
	x <- as.matrix(data)
	confusion <- NULL
	if (length(N)==0) {
		N <- 1
		p <- length(data)
		x <- matrix(data,N,p)
	}
	if (length(object$scaling)!=0){
		x <- scale(x,center=object$scaling$mu,scale=object$scaling$sd)
	}
	
	if(length(b)==1) b <- rep(b,length=K)
	if (length(a)==1) a <- matrix(a,K,max(d))
	else if (length(a)==K) a <- matrix(a,K,max(d))
	else if (object$model=='AJBQD') a <- matrix(a,K,d[1],byrow=TRUE)
	
	if(min(a,na.rm=TRUE)<=0 | min(b)<=0) stop("Some parameters A or B are negative. Prediction can't be done.\nThe reduction of the intrinsic dimensions or a more constrained model can be a solution.\nAlso you can change the value of A's and B's manually by accessing the paramaters.\n",call.=FALSE)
	
	
	if(object$model=="AJBQD") {
		K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if (object$model=="ABQD") {
		K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
	}
	else{
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- x-matrix(mu[i,],N,p,byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])
		}
	}
	
	t <- matrix(0,N,K,dimnames=list(1:N,1:K))
	for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
	result <- max.col(t)
	
	if (!is.null(object$kname)){
		result <- factor(result,labels=object$kname,levels=seq_along(object$kname))
		colnames(t) <- object$kname
	}
	if (!is.null(cls)){
		if(is.null(object$kname)) confusion <- pck_hdda_tclass(cls,result)
		else confusion <- pck_hdda_tclass(cls,result,TRUE,object$kname)
	}
	class(t) <- 'hd'
	list(class=result,posterior=t,confusion=confusion)
}

