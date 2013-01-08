hddc  <- 
function(data,K=1:10,model=c("AkjBkQkDk"),threshold=0.2,com_dim=NULL,itermax=60,eps=1e-3,graph=FALSE,algo='EM',d="Cattell",init='kmeans',show=TRUE,mini.nb=c(5,10),scaling=FALSE,min.individuals=2,noise.ctrl=1e-8,...){
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD")
	Mod2 <- c("AKJBKQKDK","AKBKQKDK ","ABKQKDK  ","AKJBQKDK ","AKBQKDK  ","ABQKDK   ","AKJBKQKD ","AKBKQKD  ","ABKQKD   ","AKJBQKD  ","AKBQKD   ","ABQKD    ","AJBQD    ","ABQD     ")
	Alg <- c('EM','CEM','SEM')
	Init <- c('random','kmeans','mini-em','param')
	algo=toupper(algo)
	if(length(init)>1){
		init <- unclass(init)
		if(any(K!=max(init))) stop("The number of class of K and of the initialization vector are different\n")
		if(length(init)!=nrow(data)) stop("The size of the initialization vector is different of the size of the data\n")
	}
	
	if(!is.numeric(d)) d <- toupper(d)
	if (d%in%c("CATTELL","C")) d <- "C"
	else if (d%in%c("BIC","B")) d <- "B"
	save_d <- d
	
	if(length(model)==1 && toupper(model)=="ALL") model <- 1:14
	if (is.numeric(model)) model <- na.omit(Mod[model])
	else model <- toupper(model)
	num <- which(model%in%as.character(1:14))
	model[num] <- Mod[as.numeric(model[num])]
	mod_num <- c()
	for(i in 1:length(model)) mod_num[i] <- which(model[i]==Mod)
	mod_num <- sort(unique(mod_num))
	model <- Mod[mod_num]
	
	if(is.integer(d)) stop("d must be equal to \"BIC\" or \"Cattell\"")
	if (any(!model%in%Mod)) stop("Invalid model name\n")
	if (any(model%in%Mod[7:14]) && is.numeric(d) && d>ncol(data)) stop("d must be strictly inferior to the dimension, \nwhich is in this case ",ncol(data),'\n')
	if (!is.numeric(min.individuals) || min.individuals<2) stop("The minimum population control variable must be superior or equal to 2.\n")
	if (length(init)==1 && !any(init==Init)) stop("Invalid initialization name\n")
	if (is.numeric(threshold)==0 || threshold<=0 || threshold>=1) stop("The parameter 'threshold' must be a double strictly within ]0,1[\n")
	if (!any(Alg==algo)) stop("Invalid algorithm name\n")
	if (length(init)==1 && init=='param' && nrow(data)<ncol(data)) stop("The 'param' initialization can't be done when N<p\n")
	if (any(is.na(data))) stop("NA values are not supported\n")
	if (length(init)==1 && init=='param' && library(MASS,logical.return=TRUE)==FALSE) stop("You need the library MASS to use the 'param' initialization\n") 
	if (length(init)==1 && init=='mini-em' && (length(mini.nb)!=2 | is.numeric(mini.nb)!=1)) stop("The parameter mini.nb must be a vector of length 2 with integers\n")
	if (!is.numeric(K) || min(K)<1) stop("K must be a vector of positive integers\n")

	data <- as.matrix(data)
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
	}
	else scaling <- NULL
	BIC <- c()
	p <- ncol(data)
	e <- vector(mode="list",length=length(K))
	if (show) cat('\t  Model  \t   K\t   BIC\n')
	nm <- length(model)
	ind <- 1
	for (i in (K <- floor(sort(K)))){
		if (i==1){
			e[[1]] <- pck_hddc(data,1,"AKJBKQKDK",threshold,d,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,...)
			BIC[1:nm] <- pck_hdda_bic(e[[1]],p)
			if (show) cat('\t',"ALL      ",'\t',1,'\t',BIC[1],'\n')
			ind <- nm+1
		}
		else {
			for (M in model){
				e[[ind]] <- pck_hddc(data,i,M,threshold,d,algo,itermax,eps,init,mini.nb,min.individuals,noise.ctrl,com_dim,...)
				if (length(e[[ind]])==1){
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',"STOPPED: pop<min.indiv.\n")
					BIC[ind] <- -Inf
				}
				else {
					if(M%in%Mod[13:14]) BIC[ind] <- pck_hdda_bic(e[[ind]],p,data)
					else BIC[ind] <- pck_hdda_bic(e[[ind]],p)
					if (show) cat('\t',Mod2[which(Mod==M)],'\t',i,'\t',BIC[ind],'\n')
				}
				ind <- ind+1
			}
		}
	}
	
	if (max(BIC)==-Inf) return(NULL)
	if(graph){
		g <- matrix(BIC,nm,length(K))
		g[g==-Inf] <- NA
		if (length(K)==1) plot(as.factor(model),g,ylab="BIC",main=paste("K =",K),xlab="model")
		else{
			plot(K,g[1,],type='o',ylim=c(min(g,na.rm=TRUE),max(BIC)),pch=1,ylab="BIC")
			if (nm>1) for (i in 2:nm) lines(K,g[i,],col=i,pch=i,type='o',lty=i)
			legend(min(K,na.rm=TRUE),max(BIC),model,col=1:nm,pch=1:nm,bty="n",lwd=1,cex=0.85,lty=1:nm)
		}
	}
	prms <- e[[which.max(BIC)]]
	if (show & (length(model)>1 | length(K)>1)) cat("\nSELECTED: model ",prms$model," with ",prms$K," clusters, BIC=",max(BIC),".\n",sep="")
	prms$BIC <- max(BIC)
	prms$scaling <- scaling
	prms$threshold <- if(save_d=="C") threshold else NULL
	class(prms) <- 'hdc'
	return(prms)
}

