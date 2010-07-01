pck_hdda_prms <-
function(data,cls,model,threshold,graph,dfixed,kname){
	p<-ncol(data)
	N<-nrow(data)
	DATA<-as.matrix(data)
	n<-prop<-c()
	k<-max(cls)
	for (i in 1:k) n[i]<-sum(cls==i)
	prop<-matrix(n/N,1,k,dimnames=list(c(''),"Prior probabilities of groups :"=kname))
	
	mu<-matrix(,k,p,dimnames=list("Class"=kname,"Group means :"=paste('V',1:p,sep='')))
	for (i in 1:k) mu[i,]<-colMeans(DATA[cls==i,])

	#calcul des matrices de variances/covariances et des vecteurs propres
	
	if (N<p) {
		if(model=="AIBQD" | model=="ABQD"){
			Y<-matrix(0,N,p)
			for (i in 1:k) Y[which(cls==i),]<-(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(N)
			donnees<-eigen(tcrossprod(Y),symmetric=TRUE)
			ev<-donnees$values
		}
		else {
			Y<-vector(mode='list',length=k)
			ev<-matrix(0,k,max(n))
			Q<-vector(mode='list',length=k)
			for (i in 1:k){ 
				Y[[i]]<-(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/sqrt(n[i])
				donnees<-eigen(tcrossprod(Y[[i]]),symmetric=TRUE)
				ev[i,1:n[i]]<-donnees$values
				Q[[i]]<-donnees$vectors
			}
		}
	}
	else if (model=="AIBQD" | model=="ABQD"){
		W<-matrix(0,p,p)
		for (i in 1:k) W<-W+prop[i]*crossprod(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i]
		donnees<-eigen(W,symmetric=TRUE)
		ev<-donnees$values
	}
	else{
		ev<-matrix(0,k,p)
		Q<-vector(mode='list',length=k)
		for (i in 1:k){ 
			donnees<-eigen(crossprod(DATA[which(cls==i),]-matrix(mu[i,],sum(cls==i),p,byrow=TRUE))/n[i],symmetric=TRUE)
			ev[i,]<-donnees$values
			Q[[i]]<-donnees$vectors
		}
	}
	
	#calcul des dimensions + graphiques
	
	if (model=='AKIBKQKD'|model=='AKBKQKD'|model=='ABKQKD'|model=='AKIBQKD'|model=='AKBQKD'|model=='ABQKD') {
		d<-rep(dfixed,length=k)
		if (graph==TRUE) {
			nbis<-ncol(ev)-1
			x<-abs(t(diff(t(ev[,1:min(d[i]+21,nbis)]))))
			x11()
			par(mfrow=c(k*(k<=4)+4*(k>4),2*(1*(k%%4!=0)+floor(k/4))))
			for (i in 1:k){	
				sub1<-paste("Class #",i,", d",i,"=",d[i],sep="")
				plot(ev[i,1:(min(d[i]+10,nbis))],type="h",col="green",main=paste("Ordered Eigen Values\nClass #",i,sep=""),xlab="",ylab="",lwd=3)
				plot(x[i,1:(min(d[i]+10,nbis-1))],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab="",xlab="")
				points(d[i],x[i,d[i]],col='red')
			}
		}
	}
	else if (model=="AIBQD" | model=="ABQD"){
		if (dfixed!='B'){
			if (length(dfixed)==1 && is.numeric(dfixed)==1) d<-rep(dfixed,k)
			else{
				x<-abs(diff(ev))
				Nmax<-max(min(floor(N/k-2),p-2),2)
				for (j in Nmax:1) {
					maxi<-max(x[1:Nmax])
					if (x[j]>=threshold*maxi) {
						d<-j
						break
					} 
					else if(j==1) d<-1
				}
				d<-rep(d,k)
			}
			if (graph) {
				nbis<-length(ev)-1
				x<-abs(diff(ev))
				x11()
				par(mfrow=c(2,1))
				plot(ev[1:(min(d[1]+10,nbis))],type="h",col="green",main="Ordered Eigen Values",xlab="",ylab="",lwd=3)
				plot(x[1:(min(d[1]+10,nbis-1))],type="l",col="blue",main=paste("Cattell's Scree-Test\nd=",d[1],sep=''),ylab="",xlab=paste('threshold=',threshold))
				abline(h=threshold*max(x),lty=3)	
				points(d[1],x[d[1]],col='red')
			}
		}
		else {
			d<-0
			if (graph)	x11()
			ev[ev<1e-10]=1e-10
			B<-c()
			for (kdim in 1:(min(p,n[i])-1)){
				if (d!=0 & kdim>d+10) break
				a<-sum(ev[1:kdim])/kdim
				b<-sum(ev[(kdim+1):min(N,p)])/(p-kdim)
				if (b<1e-10) b<-1e-10
				L2<--1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))
				B[kdim]<--2*L2+(p+kdim*(p-(kdim+1)/2)+1)*log(N)/N
				if ( (d==0 & kdim>1) && B[kdim-1]<B[kdim] )	d<-kdim-1
			}
			if (d==0) d<-min(n)-1
			if (graph){
				plot(-B,type='l',col=4,main=paste("BIC evolution w.r.t. the dimension\nd=",d,sep=''),ylab='BIC',xlab="Dimension")
				points(d,-B[d],col=2)
			}		
			d<-rep(d,k)	
		}
		if (d[1]>min(n)-1) d=rep(min(n)-1)
	}
	else if (dfixed=="B"){
		d<-rep(0,k)
		if (graph){
			x11()
			par(mfrow=c(k*(k<=4)+4*(k>4),1*(k%%4!=0)+floor(k/4)))
		}
		ev[ev<1e-10]=1e-10
		for (i in 1:k) {
			B<-c()
			maxi<-min(p,n[i])-1
			for (kdim in 1:maxi){
				if ((d[i]!=0 & kdim>d[i]+10)) break
				a<-sum(ev[i,1:kdim])/kdim
				b<-sum(ev[i,(kdim+1):(maxi+1)])/(p-kdim)
				if (b<1e-10) b<-1e-10
				L2<--1/2*(kdim*log(a)+(p-kdim)*log(b)-2*log(prop[i])+p*(1+1/2*log(2*pi)))
				B[kdim]<--2*L2+(p+kdim*(p-(kdim+1)/2)+1)*log(n[i])/n[i]
				if ((d[i]==0 & kdim>1) && B[kdim-1]<B[kdim] )d[i]<-kdim-1
			}
			if (d[i]==0) d[i]<-min(n[i]-1,p-1)
			if (graph){
				plot(-B,type='l',col=4,main=paste("class #",i,", d=",d[i],sep=''),ylab='BIC',xlab="Dimension")
				points(d[i],-B[d[i]],col=2)
			}
		}
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

		if (graph==TRUE) {
			x11()
			par(mfrow=c(k*(k<=4)+4*(k>4),2*(1*(k%%4!=0)+floor(k/4))))
			for (i in 1:k){
				sub1<-paste("Class #",i,", d",i,"=",d[i],sep="")
				plot(ev[i,1:(min(d[i]+10,nbis))],type="h",col="green",main=paste("Ordered Eigen Values\nClass #",i,sep=""),xlab="",ylab="",lwd=3)
				plot(x[i,1:(min(d[i]+10,Nmax[i]))],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab="",xlab="")
				abline(h=threshold*maxi[i],lty=3)	
				points(d[i],x[i,d[i]],col='red')
			}
		}
	}
	
	#mise en place des matrices Qi
	
	if (model=="AIBQD" |model=="ABQD"){
		if (N>=p) Q<-matrix(donnees$vectors[,1:d[1]],p,d[1])
		else {
			Q<-matrix(t(Y)%*%donnees$vectors[,1:d[1]],p,d[1])
			normalise<-c()
			for(i in 1:d[1]) normalise[i]<-as.double(crossprod(Q[,i]))
			Q<-Q/matrix(sqrt(normalise),p,d[1],byrow=TRUE)		
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
	
	if (model=='AKIBKQKDK' | model=='AKIBQKDK' | model=='AKIBKQKD' | model=='AKIBQKD' ){
		ai<-matrix(NA,k,max(d),dimnames=list("Class"=kname,"Aki :"=paste("a",1:max(d),sep='')))
		for (i in 1:k) ai[i,1:d[i]]<-ev[i,1:d[i]]
	}
	else if (model=='AKBKQKDK' | model=='AKBQKDK' | model=='AKBKQKD' | model=='AKBQKD'){
		ai<-matrix(NA,1,k,dimnames=list(c("Ak :"),kname))
		for (i in 1:k) ai[i]<-sum(ev[i,1:d[i]])/d[i]
	}
	else if (model=="AIBQD"){
		ai<-matrix(ev[1:d[1]],1,d[1],dimnames=list(c("Ai :"),paste('a',1:d[1],sep='')))
	}
	else if (model=="ABQD"){
		ai<-matrix(sum(ev[1:d[1]])/d[1],dimnames=list(c("A :"),c('')))
	}
	else {
		a<-eps<-0
		for (i in 1:k) {
			a<-a+sum(ev[i,1:d[i]])*prop[i]
			eps<-eps+prop[i]*d[i]
		}
		ai<-matrix(a/eps,dimnames=list(c("A :"),c('')))
	}
	
	if (model=='AKIBKQKDK'|model=='AKBKQKDK'|model=='ABKQKDK'|model=='AKIBKQKD'|model=='AKBKQKD'|model=='ABKQKD'){
		bi<-matrix(NA,1,k,dimnames=list(c("Bk :"),kname))
		if (N>=p) for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):p])/(p-d[i])
		else for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):n[i]])/(p-d[i])
	}
	else if (model=="ABQD" | model=="AIBQD"){
		if (N>=p) bi<-matrix(sum(ev[(d[1]+1):p])/(p-d[1]),dimnames=list(c("B :"),c('')))
		else bi<-matrix(sum(ev[(d[1]+1):N])/(p-d[1]),dimnames=list(c("B :"),c('')))
	}
	else{
		b<-eps<-0
		if (N>=p){
			for (i in 1:k) {
				eps<-eps+prop[i]*d[i]
				b<-b+sum(ev[i,(d[i]+1):p])*prop[i]
			}
		}
		else {
			for (i in 1:k) {
				eps<-eps+prop[i]*d[i]
				b<-b+sum(ev[i,(d[i]+1):n[i]])*prop[i]
			}
		}
		bi<-matrix(b/(p-eps),dimnames=list(c("B :"),c('')))
	}
	d<-matrix(d,1,k,dimnames=list(c('dim :'),"Intrinsic dimensions of the classes :"=kname))
	class(prop)<-class(mu)<-class(ai)<-class(bi)<-class(d)<-"hd"
	list(model=model,k=k,d=d,a=ai,b=bi,mu=mu,prop=prop,ev=ev,Q=Q,kname=kname)
}

