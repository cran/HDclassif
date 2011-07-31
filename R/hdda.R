hdda  <- 
function(data,cls,model='AkjBkQkDk',graph=FALSE,d="Cattell",threshold=0.2,com_dim=NULL,show=TRUE,scaling=FALSE,cv.dim=1:10,cv.threshold=c(.001,.005,.05,1:9*0.1),cv.vfold=10,LOO=FALSE,dim.ctrl=1e-8){
	Mod <- c("AKJBKQKDK","AKBKQKDK","ABKQKDK","AKJBQKDK","AKBQKDK","ABQKDK","AKJBKQKD","AKBKQKD","ABKQKD","AKJBQKD","AKBQKD","ABQKD","AJBQD","ABQD","ALL")
	Mod2 <- c("AKJBKQKDK","AKBKQKDK ","ABKQKDK  ","AKJBQKDK ","AKBQKDK  ","ABQKDK   ","AKJBKQKD ","AKBKQKD  ","ABKQKD   ","AKJBQKD  ","AKBQKD   ","ABQKD    ","AJBQD    ","ABQD     ")
	
	if (is.numeric(model)) model <- na.omit(Mod[model])
	else model <- toupper(model)
	num <- which(model%in%as.character(1:14))
	model[num] <- Mod[as.numeric(model[num])]
	mod_num <- c()
	for(i in 1:length(model)) mod_num[i] <- which(model[i]==Mod)
	mod_num <- sort(unique(mod_num))
	model <- Mod[mod_num]
	if (any(!model%in%Mod)) stop("Invalid model name\n",call.=FALSE)
	
	if(length(d)>1) stop("d cannot be a vector.\n",call.=FALSE)
	if(!is.numeric(d)) d <- toupper(d)
	
	if (d%in%c("CATTELL","C")) d <- "C"
	else if (d%in%c("BIC","B")) d <- "B"
	
	if(is.integer(d)) stop("d must be equal to \"BIC\" or \"Cattell\"",call.=FALSE)
	if (!is.numeric(threshold)) stop("The parameter 'threshold' must a double within 0 and 1.\n",call.=FALSE)
	if (threshold<0 | threshold>1) stop("The parameter 'threshold' must be within 0 and 1.\n",call.=FALSE)
	if (any(is.na(data))) stop("NA values are not allowed\n",call.=FALSE)
	if (nrow(data)!=length(cls)) stop ("The class vector does not fit with the data\n",call.=FALSE)
	if((length(model)>1 || (length(model)==1 && model=="ALL")) && d == "CV") stop("A specific model must be chosen for choosing the threshold/dimension with cross-validation.\n",call.=FALSE)
	if(any(model%in%Mod[c(1:6)]) && !d%in%c("C","B","CV")) stop("To use this model d must have the value: 'c', 'b' or 'cv'; it also may be an integer but only for common dimension model. See help for information.\n",call.=FALSE)
	if( any(model%in%Mod[7:14]) && !d%in%c("C","B","CV") && is.null(com_dim) ) stop("d must have the value: 'c', 'b' or 'cv'; it also may be an integer but only for common dimension model. See help for information.\n",call.=FALSE)
	
	
	save_d <- if(d=="B") "B" else "C"
	cls <- as.factor(cls)
	names <- levels(cls)
	z <- unclass(cls)
	data <- as.matrix(data)
	K <- max(z)
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data,"scaled:center"),sd=attr(data,"scaled:scale"))
	}
	else scaling <- NULL
	N <- nrow(data)
	p <- ncol(data)
	n <- table(z)
	
	if( (any(model%in%Mod[7:12]) && !is.null(com_dim) && com_dim > min(n,p)-1) )  stop("com_dim has to be lower or equal to ",min(n,p)-1, if(p>min(n)) paste(" because of the class",names[which.min(n)],"which has",min(n),"observations.\n") else paste(" because there are only",p,"dimensions.\n"),call.=FALSE)
	if(any(model%in%Mod[13:15]) && !is.null(com_dim) && com_dim > min(N,p)-1) stop("com_dim has to be lower or equal to ",min(N,p)-1, if(p<=N) paste(" because there are only",p,"dimensions.\n") else paste(" because there are only",N,"observations.\n"),call.=FALSE)
	
	
	if(LOO){
		if(d=="CV" && is.null(com_dim)) stop("To do LOO, d must be equal to 'c' or 'b', or you may select manually the dimension for common dimension models using 'com_dim'.\n",call.=FALSE)
		if(length(model)>1 || model=="ALL") stop("A specific model must be chosen for LOO.\n",call.=FALSE)
		
		if(model%in%Mod[1:14] && !is.null(com_dim)) d <- com_dim
		
		class <- rep(NA,N)
		posterior <- matrix(NA,N,K)
		for(i in 1:N){
			prms <- NULL
			try(prms <- pck_hdda_prms(data[-i,],z[-i],model,threshold,d,names,dim.ctrl),silent=TRUE)
			if(!is.null(prms)) {
				res <- NULL
				try(res <- predict.hdc(prms,data[i,]),silent=TRUE)
				if(!is.null(res)){
					class[i] <- res$class
					posterior[i,] <- res$posterior
				}
			}
		}
		class <- factor(class,labels=names,levels=seq_along(names))
		return(list(class=class,posterior=posterior))
	}
	
	
	if(d=='CV'){
		d.max <- if(model%in%Mod[7:12]) min(n,p)-1 else min(N,p)-1
		cv.dim <- sort(cv.dim,decreasing=TRUE)
		cv.dim <- cv.dim[cv.dim<=d.max]
		if(length(cv.dim)==0) stop("cv.dim must be an integer stricly inferior \nto the dimension.\n",call.=FALSE)
		cv.threshold <- sort(cv.threshold) 
		cv.threshold <- cv.threshold[cv.threshold>=0  & cv.threshold<=1]
		if(length(cv.threshold)==0) stop("cv.threshold must be a float within 0 and 1.\n",call.=FALSE)
		cv.vfold <- if(cv.vfold<N) cv.vfold else N
		
		u <- sample(1:N)
		ind <- c()
		for(i in 1:cv.vfold) ind[i] <- if(i==1) floor(N/cv.vfold) else floor((N-sum(ind))/(cv.vfold+1-i))
		fin <- cumsum(ind)
		debut <- c(1,fin[-cv.vfold]+1)
		
		if(model%in%Mod[7:14]) {
			n_cv <- length(cv.dim)
			cv.threshold <- rep(.5,n_cv)
		}
		else{
			n_cv <- length(cv.threshold)
			cv.dim <- rep("C",n_cv)
		}
		
		res <- fails <- rep(0,n_cv)
		N2 <- rep(N,n_cv)
		for(j in 1:cv.vfold){
			ind <- u[debut[j]:fin[j]]
			prms <- NULL
			i <- 0
			while((i <- i+1)<=n_cv && is.null(prms) ){
				try(prms <- pck_hdda_prms(data[-ind,],z[-ind],model,cv.threshold[i],cv.dim[i],names,dim.ctrl),silent=TRUE)
				if(!is.null(prms)) try(res[i] <-  res[i] + sum(predict.hdc(prms,data[ind,])$class==cls[ind]), silent=TRUE)
				else {
					N2[i] <- N2[i]-length(ind)
					fails[i] <- fails[i] + 1
				}
			}
			
			if(i<=n_cv) for(i in i:n_cv){
				if(model%in%Mod[1:6]) d <- pck_hdclassif_dim_choice(prms$ev,as.vector(table(z[-ind])),"C",cv.threshold[i],FALSE,dim.ctrl)
				else d <- rep(cv.dim[i],K)
				if(model%in%Mod[13:14]) prms$Q <- prms$Q[,1:d[1]]
				else for(ii in 1:K) if(prms$d[ii]>1) prms$Q[[ii]] <- prms$Q[[ii]][,1:d[ii]]
				prms$d <- d	
				prms_bis <- pck_hdda_prms_bis(model,prms,p)
				try(res[i] <-  res[i] + sum(predict.hdc(prms_bis,data[ind,])$class==cls[ind]), silent=TRUE)
			}
		}
		
		if(show){
			if(model%in%Mod[7:14]) cat("\t  Model   \t dim\t CV\n")
			else cat("\t  Model   \tthreshold\t CV\n")
			for(i in n_cv:1){
				if(model%in%Mod[7:14]) cat('\t',Mod2[model==Mod],'\t',cv.dim[i],"\t",res[i]/N2[i]*100,if(fails[i]>0) paste("  Info: failed",fails,"times"),'\n')
				else cat('\t',Mod2[model==Mod],'\t',cv.threshold[i],"\t\t",res[i]/N2[i]*100,if(fails[i]>0) paste("  Info: failed",fails,"times"),'\n')
			}
		}
		
		res <- res/N2*100
		res <- res[n_cv:1]
		cv.dim <- cv.dim[n_cv:1]
		cv.threshold <- cv.threshold[n_cv:1]
		if(model%in%Mod[7:14]){
			d <- com_dim <- cv.dim[which.max(res)]
			if(show) cat("Best dimension with respect to the CV results: ",d,".\n",sep="")
			if(graph){
				barplot(res-100/K,names.arg=cv.dim,offset=100/K,col="blue", xlab="Dimensions",ylab="Correct classification rate",axes=FALSE, main=paste("Cross-Validation\n(chosen dim=",d,")",sep=""))
				axis(2,at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
		}
		else{
			d <- "C"
			threshold <- cv.threshold[which.max(res)]
			if(show) cat("Best threshold with respect to the CV results: ",threshold,".\n",sep="")
			if(graph){
				barplot(res-100/K,names.arg=cv.threshold,offset=100/K,col="blue", xlab="Thresholds",ylab="Correct classification rate",axes=FALSE, main=paste("Cross-Validation\nthreshold=",threshold,sep=""))
				axis(2,at=floor(100/K+(max(res)-100/K)/5*0:5))
			}
		}
	}	
	
	if(length(model)>1){
		e <- vector(mode="list",length=max(mod_num))
		BIC <- c()
		
		for(i in mod_num){
			e[[i]] <- pck_hdda_prms(data,z,Mod[i],threshold,d,names,dim.ctrl,com_dim)
			BIC[i] <- pck_hdda_bic(e[[i]],p)
		}
		
		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC,na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		if(show){
  			cat(" # :\t  Model  \t     BIC\n")
   			for(i in mod_num){
   				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6,na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i,':\t',Mod2[i],'\t',BIC[i],wng,'\n')
   			}
			cat("\nSELECTED: Model ",prms$model,", BIC=",prms$BIC,".\n",sep="")
   		}
		
		if(graph){			
			BIC <- BIC[!is.na(BIC)]
			min_b=min(BIC[BIC!=-Inf])
			max_b=max(BIC,na.rm=TRUE)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b,names.arg=mod_num,offset=min_b,col="blue", xlab="models",ylab="BIC",axes=FALSE, main=paste("BIC for all models\n(chosen model=",prms$model,")",sep=""))
			axis(2,at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		class(prms) <- 'hdc'
		return(prms)
	}
	else if(model=="ALL"){
		e <- vector(mode="list",length=14)
		BIC <- c()
		
		#models with var dim
		e[[1]] <- pck_hdda_prms(data,z,Mod[1],threshold,d,names,dim.ctrl)
		for (i in 2:6) e[[i]] <- pck_hdda_prms_bis(Mod[i],e[[1]],p)

		#models with common dim	
		e[[7]] <- pck_hdda_prms(data,z,Mod[7],threshold,d,names,dim.ctrl,com_dim)
		for (i in 8:12) e[[i]] <- pck_hdda_prms_bis(Mod[i],e[[7]],p)
		
		#models 13 and 14: common var/covar matrix
		e[[13]] <- pck_hdda_prms(data,z,Mod[13],threshold,d,names,dim.ctrl,com_dim)
		e[[14]] <- pck_hdda_prms_bis(Mod[14],e[[13]],p)
		
		#BIC calculation
		for(i in 1:14) BIC[i] <- pck_hdda_bic(e[[i]],p)
   
   		prms <- e[[which.max(BIC)]]
		prms$BIC <- max(BIC,na.rm=TRUE)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		if(show){
  			cat(" # :\t  Model  \t     BIC\n")
   			for(i in 1:14){
   				if(i<10) cat(' ')
				wng <- if(any(e[[i]]$b<10e-6) | any(e[[i]]$a<10e-6,na.rm=TRUE)) "info: b < 10e-6" else ""
				cat(i,':\t',Mod2[i],'\t',BIC[i],wng,'\n')
   			}
			cat("\nSELECTED: Model ",prms$model,", BIC=",prms$BIC,".\n",sep="")
   		}
		if(graph){
			min_b <- min(BIC[BIC!=-Inf])
			max_b <- max(BIC)
			BIC[BIC==-Inf] <- min_b
			barplot(BIC-min_b,names.arg=1:14,offset=min_b,col="blue", xlab="models",ylab="BIC",axes=FALSE, main=paste("BIC for all models\n(chosen model=",prms$model,")",sep=""))
			axis(2,at=floor(min_b+(max_b-min_b)/5*0:5))
		}
		class(prms) <- 'hdc'
		return(prms)
	}
	else {
		prms <- pck_hdda_prms(data,z,model,threshold,d,names,dim.ctrl,com_dim)
		prms$BIC <- pck_hdda_bic(prms,p)
		prms$scaling <- scaling
		prms$threshold <- if(save_d=="C") threshold else NULL
		
		class(prms) <- 'hdc'
		return(prms)
	}
}