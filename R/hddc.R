hddc <-
function(data,k=1:10,model='AkjBkQkDk',threshold=0.2,itermax=60,eps=1e-2,graph=FALSE,algo='EM',d=NULL,init='kmeans',show=TRUE,mini.nb=c(5,10),scaling=FALSE,ctrl=1,...){
	Mod<-c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	Mod2<-c("AKJBKQKDK","AKBKQKDK ","ABKQKDK  ","AKJBQKDK ","AKBQKDK  ","ABQKDK   ","AKJBKQKD ","AKBKQKD  ","ABKQKD   ","AKJBQKD  ","AKBQKD   ","ABQKD    ","AJBQD    ","ABQD     ")
	Alg<-c('EM','CEM','SEM')
	Init<-c('random','kmeans','mini-em','param')
	algo=toupper(algo)
	model=toupper(model)
	for (i in 1:length(model)) {
		if (!any(Mod==model[i])) stop("Error : invalid model name\n")
		if (any(model[i]==Mod[7:14]) && !is.null(d) && d>ncol(data)) stop("Error: d must be strictly inferior to the dimension, \nwhich is in this case ",ncol(data)-1,'\n')
	}
	if (!is.numeric(ctrl) || ctrl<0) cat("Error : the control variable must be a strictly positive double\n")
	else if (!any(init==Init)) cat("Error : invalid initialisation name\n")
	else if (is.numeric(threshold)==0 || threshold<=0 || threshold>=1) cat("Error : the parameter 'threshold' must be a double strictly within ]0,1[\n")
	else if (!any(Alg==algo)) cat("Error : invalid algorithm name\n")
	else if (init=='param' & nrow(data)<ncol(data)) cat("The 'param' initialisation can't be done when N<p\n")
	else if (any(is.na(data))) cat("Error : NA values are not supported\n")
	else if (init=='param' && library(MASS,logical.return=TRUE)==FALSE) cat("You need the library MASS to use the 'param' initialisation\n") 
	else if (init=='mini-em' && (length(mini.nb)!=2 | is.numeric(mini.nb)!=1)) cat("Error : the parameter mini.nb must be a vector of length 2 with integers\n")
	else if (typeof(init)!="character" && length(init)!=nrow(data)) cat("Error : length of the class must fit the data\n")
	else if (length(k)>20) cat("Error : more than 20 different classes can't be tested\n")
	else if (is.numeric(k)) {
		data<-as.matrix(data)
		if (scaling) {
			data<-scale(data)
			scaling<-list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
		}
		else scaling<-NULL
		BIC<-c()
		e<-vector(mode="list",length=length(k))
		if (show) cat('\t  Model  \t k\t   BIC\n')
		ind<-1
		for (i in sort(k)){
			for (M in model){
				e[[ind]]<-pck_hddc(data,i,M,threshold,d,graph,algo,itermax,eps,init,mini.nb,ctrl,...)
				if (length(e[[ind]])==1){
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',"STOPPED\n")
					BIC[ind]<--Inf
				}
				else {
					BIC[ind]<-pck_hdda_bic(data,M,e[[ind]])
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',BIC[ind],'\n')
				}
				ind<-ind+1
			}
		}
		if (max(BIC)==-Inf) return(NULL)
		prms<-e[[which.max(BIC)]]
		if (show & (length(model)>1 | length(k)>1)) cat("\nSELECTED : model ",prms$model," with ",prms$k," clusters.\n",sep="")
		prms$BIC<-max(BIC)
		prms$scaling<-scaling
		class(prms)<-'hdc'
		prms
	}
}

