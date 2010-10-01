pck_hdda_prms_bis <-
function(data,model,par){
	N<-nrow(data)
	p<-ncol(data)
	k<-par$k
	ev<-par$ev
	d<-par$d
	kname<-par$kname
	prop<-par$prop
	n<-prop*N
	
	if (model=='AKJBKQKDK' | model=='AKJBQKDK' | model=='AKJBKQKD' | model=='AKJBQKD' ){
		ai<-matrix(NA,k,max(d),dimnames=list("Class"=kname,"Akj :"=paste("a",1:max(d),sep='')))
		for (i in 1:k) ai[i,1:d[i]]<-ev[i,1:d[i]]
	}
	else if (model=='AKBKQKDK' | model=='AKBQKDK' | model=='AKBKQKD' | model=='AKBQKD'){
		ai<-matrix(NA,1,k,dimnames=list(c("Ak :"),kname))
		for (i in 1:k) ai[i]<-sum(ev[i,1:d[i]])/d[i]
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
	
	if (model=='AKJBKQKDK'|model=='AKBKQKDK'|model=='ABKQKDK'|model=='AKJBKQKD'|model=='AKBKQKD'|model=='ABKQKD'){
		bi<-matrix(NA,1,k,dimnames=list(c("Bk :"),kname))
		if (N>=p) for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):p])/(p-d[i])
		else for(i in 1:k) bi[i]<-sum(ev[i,(d[i]+1):n[i]])/(p-d[i])
	}
	else if (model=="ABQD" | model=="AJBQD"){
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
	class(ai)<-class(bi)<-"hd"

	list(model=model,k=k,d=d,a=ai,b=bi,mu=par$mu,prop=par$prop,ev=ev,Q=par$Q,kname=par$kname)
}

