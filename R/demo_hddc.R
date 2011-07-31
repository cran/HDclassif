demo_hddc  <- 
function(){
	data(Crabs)
	devAskNewPage(ask = FALSE)
	X <- as.matrix(Crabs[,-1])
	clx <- Crabs[,1]
	algorithm <- c("EM","CEM","SEM")
	initialization <- c("kmeans","random","param")
	model <- c("AKBKQKDK","ABQKD","AJBQD")
	while (1){
		while(!any((algo <- readline("Choose the algorithm:\n1: EM; 2: CEM; 3: SEM; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(algo)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((init <- readline("Choose initialization:\n1: kmeans; 2: random; 3: param; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(init)%in%c("q","quit"))) return(invisible())
		}
		
		while(!any((mod <- readline("Choose the model:\n1: AkBkQkDk; 2: ABQkD; 3: AjBQD; q or quit: exit \n"))==as.character(1:3))){
			if(any(tolower(mod)%in%c("q","quit"))) return(invisible())
		}
		cat("hddc(data,classes,model=\"",model[as.numeric(mod)],"\",init=\"",initialization[as.numeric(init)],"\",algo=\"",algorithm[as.numeric(algo)],"\")\n",sep="")
		demo_hddc_crabs(X,4,init=initialization[as.numeric(init)],algo=algorithm[as.numeric(algo)],model=model[as.numeric(mod)],ZZ=clx)
	}
}

