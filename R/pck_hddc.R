pck_hddc <-
function(DATA,k,model,threshold,dfixed,graph,algo,itermax,eps,init,mini.nb,ctrl,...){ 
	Mod<-c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	p<-ncol(DATA)
	N<-nrow(DATA)
	if ( any(model==Mod[7:14]) && length(dfixed)==0 ) dfixed<-pck_hddc_dim(DATA,k,threshold)
	if (k>1){
		t<-matrix(0,N,k)
		if(is.numeric(init)==1 | length(init)>1) {
			name<-unique(init)
			for (i in 1:k) t[which(init==name[i]),i]<-1
		}
		else if (init=='param'){
			MU<-colMeans(DATA)
			prop<-rep(1/k,k)
			S<-crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N
			donnees<-eigen(S,symmetric=TRUE)
			ev<-donnees$values
			if (length(dfixed)==0 && !any(model==Mod[7:14])){
				x<-abs(diff(ev))
				Nmax<-max(min(floor(N/k-2),p-2),2)
				for (j in Nmax:1) {
					if (x[j]>=threshold*max(x[1:Nmax])) {
						d<-j
						break
					} 
					else if(j==1) d<-1
				}
			}
			else d<-dfixed
			a<-ev[1:d]
			b<-sum(ev[(d[1]+1):p])/(p-d[1])
			
			Q<-donnees$vectors[,1:d]
			mu<-mvrnorm(k,MU,S)
			
			K<-diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
			
			t<-matrix(0,N,k)
			for (i in 1:k) t[,i]=1/rowSums(exp((K[i,]-t(K))/2))
		}
		else if (toupper(init)=='KMEANS') {
			mc<-match.call(expand.dots = FALSE)$...
			if (is.null(mc$algorithm)) alg="Hartigan-Wong"
			else alg=mc$algorithm
			if (is.null(mc$iter.max)) im=100
			else im=mc$iter.max
			if (is.null(mc$nstart)) nst=4
			else nst=mc$nstart
			cluster<-kmeans(DATA,k,iter.max=im,nstart=nst,algorithm=alg)$cluster
			for (i in 1:k) t[which(cluster==i),i]<-1
		}
		else if (init=='mini-em'){
			prms_best<-1
			for (i in 1:mini.nb[1]){
				prms<-pck_hddc(DATA,k,model,threshold,dfixed,FALSE,algo,mini.nb[2],eps=0,init='random',mini.nb,ctrl)
				if(length(prms)!=1){
					if (length(prms_best)==1) prms_best<-prms
					else if (prms_best$loglik[length(prms_best$loglik)]<prms$loglik[length(prms$loglik)]) prms_best<-prms
				}
			}
			if (length(prms_best)==1) return(1)
			t<-prms_best$posterior
		}
		else t<-t(rmultinom(N,1,rep(1/k,k)))
	}
	else t<-matrix(1,N,1)
	
	likely<-c()
	I<-0
	test<-Inf
	while (I<itermax && test>eps){
		I<-I+1
		if (algo!='EM' && I!=1) t<-t2
		if (k>1 && (any(is.na(t)) || any(colSums(t>1/k)<=ctrl*N/100))) return(1)
		m<-pck_hddc_m_step(DATA,k,t,model,threshold,dfixed)
		t<-pck_hddc_e_step(DATA,m)
		L<-t$L
		t<-t$t
		if (algo=='CEM') {
			t2<-matrix(0,N,k)
			t2[cbind(1:N,max.col(t))]<-1
		}
		else if(algo=='SEM') { 
			t2<-matrix(0,N,k)
			for (i in 1:N)	t2[i,]<-t(rmultinom(1,1,t[i,]))
		}
		likely[I]<-L
		if (I!=1) test<-abs(likely[I]-likely[I-1])
	}

	if (graph==TRUE) {
		x11()
		plot(likely,type='l',col=2,ylab='log Likelihood',main="Log Likelihood Evolution",xlab=paste("Model = ",model," ; Algorithm =",algo,sep=""))
	}
	
	if (model=='AKBKQKDK' | model=='AKBQKDK' | model=='AKBKQKD' | model=='AKBQKD') {
		a<-matrix(m$a[,1],1,m$k,dimnames=list(c("Ak :"),1:m$k))
	}
	else if(model=='AJBQD') {
		a<-matrix(m$a[1,],1,m$d[1],dimnames=list(c('Aj :'),paste('a',1:m$d[1],sep='')))
	}
	else if (model=='ABKQKDK' | model=='ABQKDK' | model=='ABKQKD' | model=='ABQKD'|model=="ABQD") {
		a<-matrix(m$a[1],dimnames=list(c('A :'),c('')))
	}
	else a<-matrix(m$a,m$k,max(m$d),dimnames=list('Class'=1:m$k,paste('a',1:max(m$d),sep='')))
	
	if (model=='AKJBQKDK'|model=='AKBQKDK'|model=='ABQKDK'|model=='AKJBQKD'|model=='AKBQKD'|model=='ABQKD' |model=='AJBQD'|model=="ABQD") {
		b<-matrix(m$b[1],dimnames=list(c('B :'),c('')))
	}
	else b<-matrix(m$b,1,m$k,dimnames=list(c("Bk :"),1:m$k))
	
	d<-matrix(m$d,1,m$k,dimnames=list(c('dim :'),"Intrinsic dimensions of the classes :"=1:m$k))
	mu<-matrix(m$mu,m$k,p,dimnames=list('Class'=1:m$k,'Posterior group means :'=paste('V',1:p,sep='')))
	prop<-matrix(m$prop,1,m$k,dimnames=list(c(''),'Posterior probabilities of groups'=1:m$k))
	class(b)<-class(a)<-class(d)<-class(prop)<-class(mu)<-class(t)<-'hd'
	cls<-c()
	for (i in 1:N) cls[i]<-which.max(t[i,])
	list(model=model,k=k,d=d,a=a,b=b,mu=mu,prop=prop,ev=m$ev,Q=m$Q,loglik=likely,posterior=t,class=cls)
}

