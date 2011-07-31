print.hdc <-
function(x,...){
	if(length(x$kname)!=0) cat ("HIGH DIMENSIONAL DISCRIMINANT ANALYSIS\nMODEL: ",x$model,"\n",sep='')
	else cat ("HIGH DIMENSIONAL DATA CLUSTERING\nMODEL: ",x$model,"\n",sep='')
	print(x$prop)
	print(x$d)
	print(x$a)
	print(x$b)
	cat("BIC: ",x$BIC,"\n")
	if(!is.null(x$info)) cat(x$info,"\n")
	if(min(x$a,na.rm=TRUE)<0) cat("Information: a < 0\n")
	if(min(x$b)<10e-6) cat("Information: b < 10e-6\n")
}

