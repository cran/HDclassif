pck_hdda_prms_bis  <- 
function(model,par,p){
	N <- par$N
	K <- par$K
	ev <- par$ev
	d <- par$d
	kname <- par$kname
	prop <- par$prop
	n <- prop*N
	
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
	class(ai) <- class(bi) <- "hd"

	list(model=model,K=K,d=d,a=ai,b=bi,mu=par$mu,prop=par$prop,ev=ev,Q=par$Q,kname=par$kname,info=par$info,N=N,com_ev=par$com_ev)
}

