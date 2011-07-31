pck_hdda_prms  <- 
function(data,cls,model,threshold,method,kname,dim.ctrl,com_dim=NULL){
	p <- ncol(data)
	N <- nrow(data)
	K <- max(cls)
	com_ev <- NULL
	info <- NULL
	n <- as.vector(table(cls))
	prop <- matrix(n/N,1,K,dimnames=list(c(''),"Prior probabilities of groups:"=kname))
	
	mu <- matrix(rowsum(data,cls)/n,K,p,dimnames=list("Class"=kname,"Group means:"=paste('V',1:p,sep='')))

	#calcul des matrices de variance/covariance et des vecteurs propres
	
	if( model%in%c("AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD") ){
		if (N<p) {
			Y <- matrix(0,N,p)
			for (i in 1:K) Y[which(cls==i),] <- (data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(N)
			if(model%in%c("AJBQD","ABQD")) donnees <- eigen(tcrossprod(Y),symmetric=TRUE)
			else donnees <- eigen(tcrossprod(Y),symmetric=TRUE,only.values=TRUE)
		}
		else{
			W <- matrix(0,p,p)
			for (i in 1:K) W <- W + prop[i]*crossprod(data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i]
			if(model%in%c("AJBQD","ABQD")) donnees <- eigen(W,symmetric=TRUE)
			else donnees <- eigen(W,symmetric=TRUE,only.values=TRUE)
		}	
		ev <- com_ev <- donnees$values
	}
	
	if(!model%in%c("AJBQD","ABQD")){
		if(any(n<p)) Y <- vector(mode='list',length=K) #je cree le vecteur entier, c'est plus pratique
		Q <- vector(mode='list',length=K)
		ev <- matrix(NA,K,min(max(n),p))
		for(i in which(n<p)){
			Y[[i]] <- (data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(n[i])
			donnees <- eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
			ev[i,1:n[i]] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
		for(i in which(n>=p)){
			donnees <- eigen(crossprod(data[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i],symmetric=TRUE)
			ev[i,] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}
	
	#calcul des dimensions + graphiques
	
	if(model%in%c("AJBQD","ABQD")){
		if(!is.null(com_dim)) method <- com_dim
		if(method%in%c("C","B")) method <- pck_hdclassif_dim_choice(com_ev,n,method,threshold,FALSE,dim.ctrl)
		d <- rep(method,K)
	}
	else if (model%in%c('AKJBKQKD','AKBKQKD','ABKQKD','AKJBQKD','AKBQKD','ABQKD')){
		if(!is.null(com_dim)) method <- com_dim
		if(method%in%c("C","B")) method <- pck_hdclassif_dim_choice(com_ev,n,method,threshold,FALSE,dim.ctrl)
		d <- rep(method,K)
		if( d[1]>min(n,p)-1 ) {
			d[] <- min(n,p)-1
			info <- paste("Information: d has been lowered to",d[1],"because of the class",kname[which.min(n)],"which has",min(n),"observations.")
		}
		dmax <- if(any(ev<dim.ctrl,na.rm=TRUE)) min(unlist(apply(ev<dim.ctrl,1,which)))-2 else Inf
		if(d[1] > dmax) d[] <- dmax
	}
	else d <- pck_hdclassif_dim_choice(ev,n,method,threshold,FALSE,dim.ctrl)
	
	#mise en place des matrices Qi
	
	if (model%in%c("AJBQD","ABQD")){
		if (N>=p) {
			Q <- matrix(donnees$vectors[,1:d[1]],p,d[1])
		}
		else {
			Q <- matrix(t(Y)%*%donnees$vectors[,1:d[1]],p,d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[,i]))
			Q <- Q/matrix(sqrt(normalise),p,d[1],byrow=TRUE)
		}
	}
	else{
		for(i in which(n>=p)){
			Q[[i]] <- matrix(Q[[i]][,1:d[i]],p,d[i])
		}
		for(i in which(n<p)){
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][,1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][,j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise),p,d[i],byrow=TRUE)
		}
	}
	
	#calcul des paramètres
	
	if ( model%in%c('AKJBKQKDK','AKJBQKDK','AKJBKQKD','AKJBQKD') ){
		ai <- matrix(NA,K,max(d),dimnames=list("Class"=kname,"Akj:"=paste("a",1:max(d),sep='')))
		for (i in 1:K) ai[i,1:d[i]] <- ev[i,1:d[i]]
	}
	else if ( model%in%c('AKBKQKDK','AKBQKDK' ,'AKBKQKD','AKBQKD') ){
		ai <- matrix(NA,1,K,dimnames=list(c("Ak:"),kname))
		for (i in 1:K) ai[i] <- sum(ev[i,1:d[i]])/d[i]
	}
	else if (model=="AJBQD"){
		ai <- matrix(ev[1:d[1]],1,d[1],dimnames=list(c("Aj:"),paste('a',1:d[1],sep='')))
	}
	else if (model=="ABQD"){
		ai <- matrix(sum(ev[1:d[1]])/d[1],dimnames=list(c("A:"),c('')))
	}
	else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i,1:d[i]])*prop[i]
		ai <- matrix(a/eps,dimnames=list(c("A:"),c('')))
	}
	
	if ( model%in%c('AKJBKQKDK','AKBKQKDK','ABKQKDK','AKJBKQKD','AKBKQKD','ABKQKD') ){
		bi <- matrix(NA,1,K,dimnames=list(c("Bk:"),kname))
		for(i in which(n>=p)) bi[i] <- sum(ev[i,(d[i]+1):p])/(p-d[i])
		for(i in which(n<p)) bi[i] <- sum(ev[i,(d[i]+1):n[i]])/(p-d[i])
	}
	else if ( model%in%c("ABQD","AJBQD") ){
		if (N>=p) bi <- matrix(sum(ev[(d[1]+1):p])/(p-d[1]),dimnames=list(c("B:"),c('')))
		else bi <- matrix(sum(ev[(d[1]+1):N])/(N-d[1]),dimnames=list(c("B:"),c('')))
	}
	else{
		b <- 0
		eps <- sum(prop*d)
		for(i in which(n>=p)) b <- b + sum(ev[i,(d[i]+1):p])*prop[i]
		for(i in which(n<p)) b <- b + sum(ev[i,(d[i]+1):n[i]])*prop[i]
		bi <- matrix(b/(min(max(n),p)-eps),dimnames=list(c("B:"),c('')))
	}
	d <- matrix(d,1,K,dimnames=list(c('dim:'),"Intrinsic dimensions of the classes:"=kname))
	class(prop) <- class(mu) <- class(ai) <- class(bi) <- class(d) <- class(ev) <- "hd"
	list(model=model,K=K,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q,kname=kname,info=info,N=N,com_ev=com_ev)
}

