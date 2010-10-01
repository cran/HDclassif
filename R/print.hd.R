print.hd <-
function(x,...){
	class(x)<-NULL
	print(x,digits=3,na.print='.')
	class(x)<-'hd'
}

