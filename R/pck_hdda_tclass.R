pck_hdda_tclass  <- 
function(cls,cls2,confusion=TRUE,nom=NULL){
	if(!is.null(nom)){
		u <- u2 <- as.matrix(sort(unique(nom)))
		k <- length(u)
	}
	else{
		u <- as.matrix(sort(unique(cls)))
		u2 <- 1:max(cls2)
		lu <- length(u)
		lu2 <- length(u2)
		k <- max(lu,lu2)
		if(length(u)<k) u <- c(u,paste("VAR",(lu+1):k,sep=''))
		else if (length(u2)<k) u2 <- 1:k
	}
	N <- length(cls)
	x <- y <- c()
	for (i in 1:k) {
		x[which(cls==u[i])] <- i
		y[which(cls2==u2[i])] <- i
	}
	name <- u2
	comp <- matrix(,k,k)

	xclass <- yclass <- matrix(0,N,k)
	for (i in 1:k){
		xclass[which(x==i),i] <- 1
		yclass[which(y==i),i] <- 1
	}
	for (i in 1:k) for (j in 1:k) comp[i,j] <- sum(xclass[,j]*yclass[,i])
	
	if(!is.null(nom)) total <- sum(diag(comp))
	else{
		comp2 <- comp
		total <- 0
		for (i in 1:k) {
			max_comp <- c()
			for (j in 1:k) max_comp[j] <- max(comp2[j,])
			total <- total+max(max_comp)
			A <- which.max(max_comp)
			B <- which.max(comp2[A,])
			comp2[A,] <- comp2[,B] <- -1
			if (A!=B){
				svg <- comp[A,]
				comp[A,] <- comp[B,]
				comp[B,] <- svg
				Sname <- name[A]
				name[A] <- name[B]
				name[B] <- Sname
				svg2 <- comp2[A,]
				comp2[A,] <- comp2[B,]
				comp2[B,] <- svg2
			}
		}
	}
	comp <- matrix(comp,k,k,dimnames=list('Predicted class'=name,'Initial class'=u[1:k]))
	cat("Correct classification rate: ",total/length(y),".\n",sep="")
	if (confusion){
		print(comp)
		return(comp)
	}
}

