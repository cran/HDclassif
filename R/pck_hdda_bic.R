pck_hdda_bic  <- 
function(par,p,data=NULL){
	model <- par$model
	K <- par$K
	d <- par$d
	b <- par$b
	a <- par$a
	mu <- par$mu
	N <- par$N
	prop <- par$prop
	
	if(length(b)==1){
		#mise a jour des b pour comparaison avec modeles dimension variable
		eps <- sum(prop*d)
		n_max <- if(model%in%c("ABQD","AJBQD")) length(par$ev) else ncol(par$ev)
		b <- b*(n_max-eps)/(p-eps)
		b <- rep(b,length=K)
	}	
	if (length(a)==1) a <- matrix(a,K,max(d))
	else if (length(a)==K) a <- matrix(a,K,max(d))
	else if (model=='AJBQD') a <- matrix(a,K,d[1],byrow=TRUE)
	
	if(min(a,na.rm=TRUE)<=0 | any(b<0)) return(-Inf)
	
	if (is.null(par$loglik)){
		som_a <- c()
		for (i in 1:K) som_a[i] <- sum(log(a[i,1:d[i]]))
		L <-  -1/2*sum(prop * (som_a + (p-d)*log(b) - 2*log(prop) + p*(1+log(2*pi))))*N
	}
	else if (model%in%c("ABQD","AJBQD")){
		Q <- rep(list(par$Q),K)
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- data-matrix(mu[i,],N,p,byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
		A <- -1/2*t(K_pen)
		L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	}
	else L <- par$loglik[length(par$loglik)]
	
	
	ro <- K*p+K-1
	tot <- sum(d*(p-(d+1)/2))
	D <- sum(d)
	d <- d[1]
	to <- d*(p-(d+1)/2)
	if (model=='AKJBKQKDK') m <- ro+tot+2*K+D
	else if (model=='AKBKQKDK') m <- ro+tot+3*K
	else if (model=='ABKQKDK') m <- ro+tot+2*K+1
	else if (model=='AKJBQKDK') m <- ro+tot+K+D+1
	else if (model=='AKBQKDK') m <- ro+tot+2*K+1
	else if (model=='ABQKDK') m <- ro+tot+K+2
	else if (model=='AKJBKQKD') m <- ro+K*(to+d+1)+1
	else if (model=='AKBKQKD') m <- ro+K*(to+2)+1
	else if (model=='ABKQKD') m <- ro+K*(to+1)+2
	else if (model=='AKJBQKD') m <- ro+K*(to+d)+2
	else if (model=='AKBQKD') m <- ro+K*(to+1)+2
	else if (model=='ABQKD') m <- ro+K*to+3
	else if (model=='AJBQD') m <- ro+to+d+2
	else if (model=='ABQD') m <- ro+to+3
	bic <- -2*L+m*log(N)
	return(-bic)
}

