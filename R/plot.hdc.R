plot.hdc  <- 
function(x,method=NULL,threshold=NULL,dim.ctrl=1e-8,...){
	method <- if(!is.null(method)) method else if(!is.null(x$threshold)) "C" else "B"
	method <- toupper(method)
	method <- if (method%in%c("CATTELL","C")) "C" else if (method%in%c("BIC","B")) "B" 
	threshold <- if(!is.null(threshold)) threshold else if(!is.null(x$threshold)) x$threshold else  0.2
	
	if(!method%in%c("C","B")) stop("Wrong method name.\n",call.=FALSE)
	
	k <- x$K
	N <- x$N
	n <- x$prop*N
	
	if(is.null(x$com_ev)) d <- pck_hdclassif_dim_choice(x$ev,n,method,threshold,TRUE,dim.ctrl)
	else d <- pck_hdclassif_dim_choice(x$com_ev,n,method,threshold,TRUE,dim.ctrl)
}

