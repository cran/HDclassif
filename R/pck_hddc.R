pck_hddc <-
function(DATA,K,model,threshold,method,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,com_dim=NULL,...){ 
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	com_ev <- NULL
	if ( any(model==Mod[7:14]) ){
		MU <- colMeans(DATA)
		if (N<p) {
			Y <- (DATA-matrix(MU,N,p,byrow=TRUE))/sqrt(N)
			YYt <- tcrossprod(Y)
			com_ev <- eigen(YYt,symmetric=TRUE,only.values=TRUE)$values
		}
		else{
			S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
			com_ev <- eigen(S,symmetric=TRUE,only.values=TRUE)$values
		}
		if(is.null(com_dim)) com_dim <- pck_hdclassif_dim_choice(com_ev,N,method,threshold,FALSE,noise.ctrl)
	}
	if (K>1){
		t <- matrix(0,N,K)
		if(is.numeric(init)==1 | length(init)>1) {
			name <- unique(init)
			for (i in 1:K) t[which(init==name[i]),i] <- 1
		}
		else if (init=='param'){
			MU <- colMeans(DATA)
			prop <- rep(1/K,K)
			S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
			donnees <- eigen(S,symmetric=TRUE)
			ev <- donnees$values
			d <- if(is.numeric(method)) method else pck_hdclassif_dim_choice(ev,N,method,threshold,FALSE,noise.ctrl)
			a <- ev[1:d]
			b <- sum(ev[(d[1]+1):p])/(p-d[1])
			
			Q <- donnees$vectors[,1:d]
			mu <- mvrnorm(K,MU,S)
			
			K_pen <- diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
			
			t <- matrix(0,N,K)
			for (i in 1:K) t[,i]=1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
		}
		else if (init=='kmeans') {
			mc <- match.call(expand.dots = FALSE)$...
			if (is.null(mc$algorithm)) alg="Hartigan-Wong"
			else alg=mc$algorithm
			if (is.null(mc$iter.max)) im=50
			else im=mc$iter.max
			if (is.null(mc$nstart)) nst=4
			else nst=mc$nstart
			cluster <- kmeans(DATA,K,iter.max=im,nstart=nst,algorithm=alg)$cluster
			for (i in 1:K) t[which(cluster==i),i] <- 1
		}
		else if (init=='mini-em'){
			prms_best <- 1
			for (i in 1:mini.nb[1]){
				prms <- pck_hddc(DATA,K,model,threshold,method,algo,mini.nb[2],0,'random',mini.nb,min.individuals,noise.ctrl,com_dim)
				if(length(prms)!=1){
					if (length(prms_best)==1) prms_best <- prms
					else if (prms_best$loglik[length(prms_best$loglik)]<prms$loglik[length(prms$loglik)]) prms_best <- prms
				}
			}
			if (length(prms_best)==1) return(1)
			t <- prms_best$posterior
		}
		else {
			t <- t(rmultinom(N,1,rep(1/K,K)))
			compteur=1
			while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N,1,rep(1/K,K)))
			if(min(colSums(t))<1) stop("Random initialization failed because of too many classes and too few observations")
		}
	}
	else t <- matrix(1,N,1)
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		if (algo!='EM' && I!=1) t <- t2
		if (K>1 && (any(is.na(t)) || any(colSums(t>1/K)<min.individuals))) return(1)
		m <- pck_hddc_m_step(DATA,K,t,model,threshold,method,noise.ctrl,com_dim)
		t <- pck_hddc_e_step(DATA,m)
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0,N,K)
			t2[cbind(1:N,max.col(t))] <- 1
		}
		else if(algo=='SEM') { 
			t2 <- matrix(0,N,K)
			for (i in 1:N)	t2[i,] <- t(rmultinom(1,1,t[i,]))
		}
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
	}
	
	if ( model%in%c('AKBKQKDK','AKBQKDK','AKBKQKD','AKBQKD') ) {
		a <- matrix(m$a[,1],1,m$K,dimnames=list(c("Ak:"),1:m$K))
	}
	else if(model=='AJBQD') {
		a <- matrix(m$a[1,],1,m$d[1],dimnames=list(c('Aj:'),paste('a',1:m$d[1],sep='')))
	}
	else if ( model%in%c('ABKQKDK','ABQKDK','ABKQKD','ABQKD',"ABQD") ) {
		a <- matrix(m$a[1],dimnames=list(c('A:'),c('')))
	}
	else a <- matrix(m$a,m$K,max(m$d),dimnames=list('Class'=1:m$K,paste('a',1:max(m$d),sep='')))
	
	if ( model%in%c('AKJBQKDK','AKBQKDK','ABQKDK','AKJBQKD','AKBQKD','ABQKD','AJBQD',"ABQD") ) {
		b <- matrix(m$b[1],dimnames=list(c('B:'),c('')))
	}
	else b <- matrix(m$b,1,m$K,dimnames=list(c("Bk:"),1:m$K))
	
	d <- matrix(m$d,1,m$K,dimnames=list(c('dim:'),"Intrinsic dimensions of the classes:"=1:m$K))
	mu <- matrix(m$mu,m$K,p,dimnames=list('Class'=1:m$K,'Posterior group means:'=paste('V',1:p,sep='')))
	prop <- matrix(m$prop,1,m$K,dimnames=list(c(''),'Posterior probabilities of groups'=1:m$K))
	class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- class(t) <- 'hd'
	cls <- max.col(t)
	list(model=model,K=K,d=d,a=a,b=b,mu=mu,prop=prop,ev=m$ev,Q=m$Q,loglik=likely,posterior=t,class=cls,com_ev=com_ev,N=N)
}

