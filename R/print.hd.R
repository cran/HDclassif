print.hd <-
function(x,...){
	class(x) <- NULL
	print.default(x,digits=3,na.print='.')
	class(x) <- 'hd'
}

