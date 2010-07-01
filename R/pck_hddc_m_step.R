pck_hddc_m_step <-
function(x,k,t,model,threshold,dfixed,graph){
	N<-nrow(x)
	p<-ncol(x)
	prop<-c()
	n<-colSums(t)
	prop<-n/N
	mu<-matrix(,k,p)
	for (i in 1:k) mu[i,]<-colSums(x*t[,i])/n[i]
	
	#calcul des matrices de variances/covariances
	
	if (N<p) {
		if(model=="AIBQD" | model=="ABQD"){
			Y<-matrix(0,N,p)
			for (i in 1:k) Y<-Y+(x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(N)*sqrt(t[,i])
			donnees<-eigen(tcrossprod(Y),symmetric=TRUE)
			ev<-donnees$values
		}
		else{
			Y<-vector(mode='list',length=k)
			ev<-matrix(0,k,N)
			Q<-vector(mode='list',length=k)
			for (i in 1:k){ 
				Y[[i]]<-(x-matrix(mu[i,],N,p,byrow=TRUE))/sqrt(n[i])*sqrt(t[,i])
				donnees<-eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
				ev[i,1:N]<-donnees$values
				Q[[i]]<-donnees$vectors
			}
		}
	}
	else if (model=="AIBQD" | model=="ABQD"){
		W<-matrix(0,p,p)
		for (i in 1:k) W<-W+crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/N
		donnees<-eigen(W,symmetric=TRUE)
		ev<-donnees$values
	}
	else {
		ev<-matrix(0,k,p)
		Q<-vector(mode='list',length=k)
		for (i in 1:k){ 
			donnees<-eigen(crossprod((x-matrix(mu[i,],N,p,byrow=TRUE))*sqrt(t[,i]))/n[i],symmetric=TRUE)
			ev[i,]<-donnees$values
			Q[[i]]<-donnees$vectors
		}
	}
	
	#détermination de la dimension
	
	if (model=='AKIBKQKD'|model=='AKBKQKD'|model=='ABKQKD'|model=='AKIBQKD'|model=='AKBQKD'|model=='ABQKD'|model=="AIBQD" | model=="ABQD") {
		d<-rep(dfixed,length=k)
	}
	else{ 
		x<-abs(t(diff(t(ev))))
		nbis<-ncol(x)
		d<-maxi<-Nmax<-c()
		for (i in 1:k) Nmax[i]<-max(min(floor(n[i]-2),p-2),2)
		for (i in 1:k) for (j in Nmax[i]:1) {
			maxi[i]<-max(x[i,1:Nmax[i]])
			if (x[i,j]>=threshold*maxi[i]) {
				d[i]<-j
				break
			} 
			else if(j==1) d[i]<-1
		}
	}
	
	#mise en place des matrices Qi
	
	if (model=="AIBQD" |model=="ABQD"){
		if (N>=p) Q<-matrix(donnees$vectors[,1:d],p,d)
		else {
			Q<-matrix(t(Y)%*%donnees$vectors[,1:d],p,d)
			normalise<-c()
			for(i in 1:d) normalise[i]<-as.double(crossprod(Q[,i]))
			Q<-Q/matrix(sqrt(normalise),p,d,byrow=TRUE)		
		}
	}
	else if (N>=p) for(i in 1:k) Q[[i]]<-matrix(Q[[i]][,1:d[i]],p,d[i])
	else{
		for (i in 1:k){ 
			Q[[i]]<-t(Y[[i]])%*%(Q[[i]][,1:d[i]])
			normalise<-c()
			for (j in 1:d[i]) normalise[j]<-as.double(crossprod(as.matrix(Q[[i]][,j])))
			Q[[i]]<-Q[[i]]/matrix(sqrt(normalise),p,d[i],byrow=TRUE)
		}
	}
	
	#calcul des paramètres	
	
	ai<-matrix(NA,k,max(d))
	if (model=='AKIBKQKDK' | model=='AKIBQKDK' | model=='AKIBKQKD' | model=='AKIBQKD' ){
		for (i in 1:k) ai[i,1:d[i]]<-ev[i,1:d[i]]
	}
	else if (model=='AKBKQKDK' | model=='AKBQKDK' | model=='AKBKQKD' | model=='AKBQKD'){
		for (i in 1:k) ai[i,]<-rep(sum(ev[i,1:d[i]])/d[i],length=max(d))
	}
	else if (model=="AIBQD") for (i in 1:k) ai[i,]<-ev[1:d]
	else if (model=="ABQD")	ai[]<-sum(ev[1:d])/d
	else {
		a<-eps<-0
		for (i in 1:k) {
			a<-a+sum(ev[i,1:d[i]])*prop[i]
			eps<-eps+prop[i]*d[i]
		}
		a<-a/eps
		ai<-matrix(a,k,max(d))
	}

	bi<-c()
	if (model=='AKIBKQKDK'|model=='AKBKQKDK'|model=='ABKQKDK'|model=='AKIBKQKD'|model=='AKBKQKD'|model=='ABKQKD'){
		if (N>=p) for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):p])/(p-d[i])
		else for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):N])/(p-d[i])
	}
	else if (model=="ABQD" | model=="AIBQD"){
		if (N>=p) bi[1:k]<-sum(ev[(d+1):p])/(p-d)
		else bi[1:k]<-sum(ev[(d+1):N])/(p-d)
		d<-rep(d,k)
	}
	else{
		b<-eps<-0
		if (N>=p){
			for (i in 1:k) {
				eps<-eps+prop[i]*d[i]
				b<-b+sum(ev[i,(d[i]+1):p])*prop[i]
			}
			bi[1:k]<-b/(p-eps)
		}
		else {
			for (i in 1:k) {
				eps<-eps+prop[i]*d[i]
				b<-b+sum(ev[i,(d[i]+1):N])*prop[i]
			}
			bi[1:k]<-b/(p-eps)
		}
	}

	list(model=model,k=k,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q)
}

