print.hdc <-
function(x,...){
	if(length(x$kname)!=0) cat ("HIGH DIMENSIONAL DISCRIMINANT ANALYSIS\n\nMODEL : ",x$model,"\n\n",sep='')
	else cat ("HIGH DIMENSIONAL DATA CLUSTERING\n\nMODEL : ",x$model,"\n\n",sep='')
	print(x$prop)
	cat('\n')
	print(x$d)
	cat('\n')
	print(x$a)
	cat('\n')
	print(x$b)
	cat("\nBIC : ",x$BIC,"\n")
}

