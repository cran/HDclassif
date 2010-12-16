demo_hddc_crabs <-
function(data,k=4,model='AKBKQKD',threshold=0.2,d=1,algo='EM',itermax=80,eps=1e-4,init='kmeans',ZZ=NULL,ctrl=1,...){ 
	p<-ncol(data)
	N<-nrow(data)
	DATA<-as.matrix(data)
	t<-matrix(0,N,k)
	if (init=='param'){
		MU<-colMeans(DATA)
		prop<-rep(1/k,k)
		if (N<p){
			Y<-(DATA-matrix(MU,N,p,byrow=TRUE))/sqrt(N)*sqrt(1/k)
			S<-tcrossprod(Y)
		}
		else S<-crossprod(DATA-matrix(MU,N,p,byrow=TRUE))/N/k
		donnees<-eigen(S,symmetric=TRUE)
		vp<-donnees$values
		a<-vp[1:d]
		if (N<p) b<-sum(vp[(d+1):N])/(p-d)
		else b<-sum(vp[(d+1):p])/(p-d)
		Q<-donnees$vectors[,1:d]
		mu<-mvrnorm(k,MU,S)
		K<-diag((mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a,d,d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		t<-matrix(,N,k)
		for (i in 1:N){
			for (j in 1:k) {
				s<-sum(exp(1/2*(K[j,i]-K[,i])))
				t[i,j]<-1/s
			}
		}
	}
	else if (init=='kmeans') {
		cluster<-kmeans(DATA,k,iter.max=100,nstart=5)$cluster
		for (i in 1:k) t[which(cluster==i),i]<-1
	}
	else t<-t(rmultinom(N,1,rep(1/k,k)))
	
	likely<-c()
	I<-0
	test<-10000
	while (I<itermax && test>eps){
		I<-I+1
		if (algo!='EM' && I!=1) t<-t2
		if (k>1 && (any(is.na(t)) || any(colSums(t>1/k)<=ctrl*N/100))) {
			cat("Algorithm stopped due to an empty class.\n")
			return(invisible())
		}
		m<-pck_hddc_m_step(DATA,k,t,model,threshold,d)
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
		
		classes<-c()
		for (i in 1:N) classes[i]=which.max(t[i,])
		for (oizrhg in 1:150000){}
		demo_hddc_acp(DATA,classes,m,xlab=paste('Iteration',I),ylab='',main="Clustering process",...)
		
		likely[I]<-L
		if (I!=1) test<-abs(likely[I]-likely[I-1])
	}

	cls<-c()
	for (i in 1:N) cls[i]<-which.max(t[i,])
	if (!is.null(ZZ)) pck_hdda_tclass(cls,ZZ,FALSE)
}

