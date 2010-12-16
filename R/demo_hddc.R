demo_hddc <-
function(){
	data(Crabs)
	devAskNewPage(ask = FALSE)
	X<-Crabs[,-1]
	clx<-Crabs[,1]
	k<-max(unclass(as.factor(clx)))
	algorithm<-c("EM","CEM","SEM")
	initialization<-c("kmeans","random","param")
	model<-c("AKBKQKDK","ABQKD","AJBQD")
	while (1){
		algo<-as.numeric(readline("Choose the algorithm :\n1 : EM ; 2 : CEM ; 3 : SEM ; else : exit\n"))
		if (is.na(algo) | !any(algo==1:3) | length(algo)!=1) return(invisible())
		else{
			init<-as.numeric(readline("Choose initialisation :\n1 : kmeans ; 2 : random ; 3 : param ; else : exit\n"))
			if (is.na(init) | !any(init==1:3) | length(init)!=1) return(invisible())
			else {
				mod<-as.numeric(readline("Choose the model :\n1 : AkBkQkDk ; 2 : ABQkD ; 3 : AjBQD ; else : exit\n"))
				if (is.na(mod) | !any(mod==1:3) | length(mod)!=1) return(invisible())
				else demo_hddc_crabs(X,k,init=initialization[init],algo=algorithm[algo],model=model[mod],ZZ=clx)
			}
		}
	}
}

