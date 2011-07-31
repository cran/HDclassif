pck_hddc_e_step  <- 
function(x,par){
	p <- ncol(x)
	N <- nrow(x)
	K <- par$K
	a <- par$a
	b <- par$b
	mu <- par$mu
	d <- par$d
	prop <- par$prop
	Q <- par$Q
	
	b[b<1e-6] <- 1e-6

	if(par$model=="AJBQD") {
		K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
	}
	else if(par$model=="ABQD") {
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
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
	}
	A <- -1/2*t(K_pen)
	L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	
	t <- matrix(0,N,K)
	for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
	list(t=t,L=L)
}

