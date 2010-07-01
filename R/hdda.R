hdda <-
function(data,cls,model='AkiBkQkDk',graph=FALSE,d="Cattell",threshold=0.2,show=TRUE,scaling=FALSE){
	Mod<-c("AKIBKQKDK","AKBKQKDK","ABKQKDK","AKIBQKDK","AKBQKDK","ABQKDK","AKIBKQKD","AKBKQKD","ABKQKD","AKIBQKD","AKBQKD","ABQKD","AIBQD","ABQD","DBEST","BEST")
	model=toupper(model)
	if (any(Mod==model)==0) stop("Error : invalid model name\n",.Call=FALSE)
	else if (is.numeric(threshold)==0) stop("Error : the parameter 'threshold' must a double strictly within ]0,1[\n",.Call=FALSE)
	else if (threshold<=0 | threshold>=1) stop("Error : the parameter 'threshold' must be strictly within ]0,1[\n",.Call=FALSE)
	else if (any(is.na(data))) stop("Error : NA values are not supported\n",.Call=FALSE)
	else if (nrow(data)!=length(cls)) stop ("Error : the class vector does  not fit with the data\n",.Call=FALSE)
	if (toupper(d)=="CATTELL") d="C"
	else if (toupper(d)=="BIC") d="B"
	cls=as.factor(cls)
	u=levels(cls)
	z<-unclass(cls)
	data<-as.matrix(data)
	if (scaling) {
		data=scale(data)
		scaling<-list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
	}
	else scaling<-NULL
	N=nrow(data)
	p=ncol(data)
	
	if (model=='BEST' | model=='DBEST') {
		if (model=='DBEST' && is.numeric(d) && (length(d)!=1 || d>((N<p)*(N-1)+(N>=p)*(p-1)))) cat("Error : in order to run this model, d must be an integer stricly inferior \nto the dimension (or the number of observations), which is in this case :",((N<p)*(N-1)+(N>=p)*(p-1)),"\n")
		else{
			Mod2<-c("AKIBKQKDK","AKBKQKDK ","ABKQKDK  ","AKIBQKDK ","AKBQKDK  ","ABQKDK   ","AKIBKQKD ","AKBKQKD  ","ABKQKD   ","AKIBQKD  ","AKBQKD   ","ABQKD    ","AIBQD    ","ABQD     ")
			e<-vector(mode="list",length=8)
			BIC<-c()
			if (model=='DBEST') {
				n<-6
				if (length(d)!=1 || is.numeric(d)==0) d<-pck_hdda_prms_dim(data,z,threshold,graph,d)
				graph=FALSE
			}
			else n<-0
			if (show) cat("\t  Model  \t    BIC\n")
			e[[1]]<-pck_hdda_prms(data,z,Mod[1+n],threshold,graph,d,u)
			BIC[1]<-pck_hdda_bic(data,Mod[1+n],e[[1]])
			if (show) cat('\t',Mod2[1+n],'\t',BIC[1],'\n')
			for (i in 2:6) {
				e[[i]]<-pck_hdda_prms_bis(data,Mod[i+n],e[[1]])
				BIC[i]<-pck_hdda_bic(data,Mod[i+n],e[[i]])
				if (show) cat('\t',Mod2[i+n],'\t',BIC[i],'\n')
			}
			if (model=='DBEST') {
				e[[7]]<-pck_hdda_prms(data,z,Mod[13],threshold,FALSE,d,u)
				BIC[[7]]<-pck_hdda_bic(data,Mod[13],e[[7]])
				if (show) cat('\t',Mod2[13],'\t',BIC[7],'\n')
				e[[8]]<-pck_hdda_prms_bis(data,Mod[14],e[[7]])
				BIC[[8]]<-pck_hdda_bic(data,Mod[14],e[[8]])
				if (show) cat('\t',Mod2[14],'\t',BIC[8],'\n')
			}
			prms<-e[[which.max(BIC)]]
			prms$BIC<-max(BIC)
			prms$scaling<-scaling
			class(prms)<-'hdc'
			prms
		}
	}
	else {
		if (any(model==Mod[7:14]) && is.numeric(d) && (length(d)!=1 || d>((N<p)*(N-1)+(N>=p)*(p-1)))) cat("Error : in order to run this model, d must be an integer stricly inferior \nto the dimension (or the numberr of observations), which is in this case :",((N<p)*(N-1)+(N>=p)*(p-1)),"\n")
		else{
			if (sum(model==Mod[7:12])==1 & (length(d)!=1 || is.numeric(d)==0)) {
				d<-pck_hdda_prms_dim(data,z,threshold,graph,d)
				graph=FALSE
			}
			prms<-pck_hdda_prms(data,z,model,threshold,graph,d,u)
			prms$BIC<-pck_hdda_bic(data,model,prms)
			prms$scaling<-scaling
			class(prms)<-'hdc'
			prms
		}
	}
}

