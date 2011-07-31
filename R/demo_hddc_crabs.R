demo_hddc_crabs <-
function(DATA,k=4,model='AKBKQKD',threshold=0.2,method='C',algo='EM',itermax=80,eps=1e-2,init='kmeans',ZZ=NULL,ctrl=1,dim.ctrl=1e-8,...){ 
	com_dim <- 1
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	t <- matrix(0,N,k)
	if(model%in%Mod[7:14]) method <- 1
	if (init=='param'){
		MU <- colMeans(DATA)
		prop <- rep(1/k,k)
		S <- crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
		donnees <- eigen(S,symmetric=TRUE)
		ev <- donnees$values
		d <- if(is.numeric(method)) method else pck_hdclassif_dim_choice(ev,N,method,threshold,FALSE,dim.ctrl)
		a <- ev[1:d]
		b <- sum(ev[(d[1]+1):p])/(p-d[1])
		
		Q <- donnees$vectors[,1:d]
		mu <- mvrnorm(k,MU,S)
		
		K <- diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		
		t <- matrix(0,N,k)
		for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
	}
	else if (init=='kmeans') {
		mc <- match.call(expand.dots = FALSE)$...
		if (is.null(mc$algorithm)) alg="Hartigan-Wong"
		else alg=mc$algorithm
		if (is.null(mc$iter.max)) im=50
		else im=mc$iter.max
		if (is.null(mc$nstart)) nst=4
		else nst=mc$nstart
		cluster <- kmeans(DATA,k,iter.max=im,nstart=nst,algorithm=alg)$cluster
		for (i in 1:k) t[which(cluster==i),i] <- 1
	}
	else {
		t <- t(rmultinom(N,1,rep(1/k,k)))
		compteur=1
		while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N,1,rep(1/k,k)))
		if(min(colSums(t))<1) stop("Random initialization failed because of too many classes and too few observations")
	}
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		if (algo!='EM' && I!=1) t <- t2
		if (k>1 && (any(is.na(t)) || any(colSums(t>1/k)<=ctrl*N/100))) return(1)
		m <- pck_hddc_m_step(DATA,k,t,model,threshold,method,dim.ctrl,1)
		t <- pck_hddc_e_step(DATA,m)
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0,N,k)
			t2[cbind(1:N,max.col(t))] <- 1
		}
		else if(algo=='SEM') { 
			t2 <- matrix(0,N,k)
			for (i in 1:N)	t2[i,] <- t(rmultinom(1,1,t[i,]))
		}
		
		classes<-c()
		for (i in 1:N) classes[i]=which.max(t[i,])
		demo_hddc_acp(DATA,classes,m,xlab=paste('Iteration',I),ylab='',main="Clustering process",...)
		Sys.sleep(0.05)
		
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
	}
	
	cls <- max.col(t)
	pck_hdda_tclass(ZZ,cls,FALSE)
}

