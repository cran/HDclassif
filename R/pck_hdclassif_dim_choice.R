pck_hdclassif_dim_choice <- 
function(ev,n,method,threshold,graph,dim.ctrl){
	N <- sum(n)
	prop <- n/N
	K <- if(is.matrix(ev)) nrow(ev) else 1
	if(is.matrix(ev) && K>1){
  		p <- ncol(ev)
		if(method=="C"){
			dev <- abs(apply(ev,1,diff))
			max_dev <- apply(dev,2,max,na.rm=TRUE)
			dev <- dev/rep(max_dev,each=p-1)
			d <- apply((dev>threshold)*(1:(p-1))*t(ev[,-1]>dim.ctrl),2,which.max)
			
			if(graph){
				par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9),1+floor(K/4)-1*(K==12)+1*(K==7)))
				for(i in 1:K){
					sub1 <- paste("Class #",i,", d",i,"=",d[i],sep="")
					Nmax <- max(which(ev[i,]>dim.ctrl))-1
					plot(dev[1:(min(d[i]+5,Nmax)),i],type="l",col="blue",main=paste("Cattell's Scree-Test\n",sub1,sep=""),ylab=paste("threshold =",threshold),xlab="Dimension",ylim=c(0,1.05))
					abline(h=threshold,lty=3)  	
					points(d[i],dev[d[i],i],col='red')
				}
			}
		}
		else if(method=="B"){
			d <- rep(0,K)
			if(graph) par(mfrow=c(K*(K<=3)+2*(K==4)+3*(K>4 && K<=9)+4*(K>9),1*(1+floor(K/4)-1*(K==12)+1*(K==7))))
			
			for (i in 1:K) {
				B <- c()
				Nmax <- max(which(ev[i,]>dim.ctrl))-1
				p2 <- sum(!is.na(ev[i,]))
				Bmax <- -Inf
				for (kdim in 1:Nmax){
					if ((d[i]!=0 & kdim>d[i]+10)) break
					a <- sum(ev[i,1:kdim])/kdim
					b <- sum(ev[i,(kdim+1):p2])/(p2-kdim)
					if (b<0 | a<0) B[kdim] <- -Inf
					else{
						L2 <- -1/2*(kdim*log(a)+(p2-kdim)*log(b)-2*log(prop[i])+p2*(1+1/2*log(2*pi))) * n[i]
						B[kdim] <- 2*L2 - (p2+kdim*(p2-(kdim+1)/2)+1) * log(n[i])
					}
					if ( B[kdim]>Bmax ){
						Bmax <- B[kdim]
						d[i] <- kdim
					}
				}
				
				if(graph){
					plot(B,type='l',col=4,main=paste("class #",i,", d=",d[i],sep=''),ylab='BIC',xlab="Dimension")
					points(d[i],B[d[i]],col=2)
				}
			}
		}
	}
  	else{
		ev <- as.vector(ev)
		p <- length(ev)
		if(method=="C"){
			dvp <- abs(diff(ev))
			Nmax <- max(which(ev>dim.ctrl))-1
			if (p==2) d <- 1
			else d <- max(which(dvp[1:Nmax]>=threshold*max(dvp[1:Nmax])))
			diff_max <- max(dvp[1:Nmax])
			
			if(graph){
				plot(dvp[1:(min(d+5,p-1))]/diff_max,type="l",col="blue",main=paste("Cattell's Scree-Test\nd=",d,sep=''),ylab=paste("threshold =",threshold,sep=' '),xlab='Dimension',ylim=c(0,1.05))
				abline(h=threshold,lty=3)	
				points(d,dvp[d]/diff_max,col='red')
			}
		}
		else if(method=="B"){
			d <- 0
			Nmax <- max(which(ev>dim.ctrl))-1
			B <- c()
			Bmax <- -Inf
			for (kdim in 1:Nmax){
				if (d!=0 && kdim>d+10) break
				a <- sum(ev[1:kdim])/kdim
				b <- sum(ev[(kdim+1):p])/(p-kdim)
				if (b<=0 | a<=0) B[kdim] <- -Inf
				else{
					L2 <- -1/2*(kdim*log(a)+(p-kdim)*log(b)+p*(1+1/2*log(2*pi)))*N
					B[kdim] <- 2*L2 - (p+kdim*(p-(kdim+1)/2)+1)*log(N)
				}
				if ( B[kdim]>Bmax ){
					Bmax <- B[kdim]
					d <- kdim
				}
			}
			
			if(graph){
				plot(B,type='l',col=4,main=paste("BIC criterion\nd=",d,sep=''),ylab='BIC',xlab="Dimension")
				points(d,B[d],col=2)
			}
		}
  	}
	return(d)
}