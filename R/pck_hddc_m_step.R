pck_hddc_m_step  <- 
function(x,K,t,model,threshold,method,dim.ctrl,com_dim){
	N <- nrow(x)
	p <- ncol(x)
	prop <- c()
	n <- colSums(t)
	prop <- n/N
	mu <- matrix(,K,p)
	for (i in 1:K) mu[i,] <- colSums(x*t[,i])/n[i]
	
	ind <- apply(t>0,2,which)
	n_bis <- c()
	for(i in 1:K) n_bis[i] <- length(ind[[i]])
	
	#calcul des matrices de variances/covariances
	
	if (N<p) {
		if( model%in%c("AJBQD","ABQD") ){
			Y <- matrix(0,N,p)
			for (i in 1:K) Y <- Y+(x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(N)*sqrt(t[,i])
			donnees <- eigen(tcrossprod(Y),symmetric=TRUE)
			ev <- donnees$values
		}
		else{
			Y <- vector(mode='list',length=K)
			ev <- matrix(0,K,N)
			Q <- vector(mode='list',length=K)
			for (i in 1:K){ 
				Y[[i]] <- (x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(n[i])*sqrt(t[,i])
				donnees <- eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
				ev[i,1:N] <- donnees$values
				Q[[i]] <- donnees$vectors
			}
		}
	}
	else if ( model%in%c("AJBQD","ABQD") ){
		W <- matrix(0,p,p)
		for (i in 1:K) W <- W + crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/N
		donnees <- eigen(W,symmetric=TRUE)
		ev <- donnees$values
	}
	else {
		ev <- matrix(0,K,p)
		Q <- vector(mode='list',length=K)
		for (i in 1:K){ 
			donnees <- eigen(crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/n[i],symmetric=TRUE)
			ev[i,] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}	
	
	#selection des dimensions
	
	if (model%in%c("AJBQD","ABQD")) d <- rep(com_dim,length=K)
	else if ( model%in%c("AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD") ){
		dmax <- min(apply((ev>dim.ctrl)*rep(1:ncol(ev),each=K),1,which.max))-1
		if(com_dim>dmax) com_dim <- dmax
		d <- rep(com_dim,length=K)
	}
	else d <- pck_hdclassif_dim_choice(ev,n,method,threshold,FALSE,dim.ctrl)
	
	#mise en place des matrices Qk	
	
	if ( model%in%c("AJBQD","ABQD") ){
		if (N>=p) Q <- matrix(donnees$vectors[,1:d[1]],p,d[1])
		else {
			Q <- matrix(t(Y)%*%donnees$vectors[,1:d[1]],p,d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[,i]))
			Q <- Q/matrix(sqrt(normalise),p,d,byrow=TRUE)
		}
	}
	else if (N>=p) for(i in 1:K) Q[[i]] <- matrix(Q[[i]][,1:d[i]],p,d[i])
	else{
		for (i in 1:K){ 
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][,1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][,j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise),p,d[i],byrow=TRUE)
		}
	}
	
	#calcul des paramètres	
	
	ai <- matrix(NA,K,max(d))
	if ( model%in%c('AKJBKQKDK','AKJBQKDK','AKJBKQKD','AKJBQKD') ){
		for (i in 1:K) ai[i,1:d[i]] <- ev[i,1:d[i]]
	}
	else if ( model%in%c('AKBKQKDK','AKBQKDK' ,'AKBKQKD','AKBQKD') ){
		for (i in 1:K) ai[i,] <- rep(sum(ev[i,1:d[i]])/d[i],length=max(d))
	}
	else if (model=="AJBQD") for (i in 1:K) ai[i,] <- ev[1:d[1]]
	else if (model=="ABQD")	ai[] <- sum(ev[1:d[1]])/d[1]
	else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i,1:d[i]])*prop[i]
		ai <- matrix(a/eps,K,max(d))
	}

	bi <- c()
	if ( model%in%c('AKJBKQKDK','AKBKQKDK','ABKQKDK','AKJBKQKD','AKBKQKD','ABKQKD') ){
		for(i in 1:K) bi[i] <- sum(ev[i,(d[i]+1):min(N,p)])/(p-d[i])
	}
	else if ( model%in%c("ABQD","AJBQD") ){
		bi[1:K] <- sum(ev[(d[1]+1):min(N,p)])/(min(N,p)-d[1])
	}
	else{		
		b <- 0
		eps <- sum(prop*d)
		for(i in 1:K) b <- b + sum(ev[i,(d[i]+1):min(N,p)])*prop[i]
		bi[1:K] <- b/(min(N,p)-eps)
	}

	list(model=model,K=K,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q)
}

