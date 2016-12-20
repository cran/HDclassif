# roxygen2::roxygenise(roclets = "rd")

#' High Dimensional Data Clustering
#'
#'HDDC is a model-based clustering method. It is based on the Gaussian Mixture Model and on the idea that the data lives in subspaces with a lower dimension than the dimension of the original space. It uses the Expectation - Maximisation algorithm to estimate the parameters of the model.
#'
#' @param data A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
#' @param K A vector of integers specifying the number of clusters for which the BIC and the parameters are to be calculated; the function keeps the parameters which maximises the \code{criterion}. Default is 1:10.
#' @param model A character string vector, or an integer vector indicating the models to be used. The available models are: "AkjBkQkDk" (default), "AkBkQkDk", "ABkQkDk", "AkjBQkDk", "AkBQkDk", "ABQkDk", "AkjBkQkD", "AkBkQkD", "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD", "AjBQD", "ABQD". It is not case sensitive and integers can be used instead of names, see details for more information. Several models can be used, if it is, only the results of the one which maximizes the BIC criterion is kept. To run all models, use model="ALL".
#' @param threshold A float stricly within 0 and 1. It is the threshold used in the Cattell's Scree-Test.
#' @param criterion Either \dQuote{BIC} or \dQuote{ICL}. If several models are run, the best model is selected using the criterion defined by \code{criterion}.
#' @param com_dim It is used only for common dimensions models. The user can give the common dimension he wants. If used, it must be an integer. Its default is set to NULL.
#' @param itermax The maximum number of iterations allowed. The default is 200.
#' @param eps A positive double. It is the stopping criterion: the algorithm stops when the difference between two successive Log Likelihoods is lower than \sQuote{eps}.
#' @param algo A character string indicating the algorithm to be used. The available algorithms are the Expectation-Maximisation ("EM"), the Classification E-M ("CEM") and the Stochastic E-M ("SEM"). The default algorithm is the "EM".
#' @param d_select Either \dQuote{Cattell} (default) or \dQuote{BIC}. See details for more information. This parameter selects which method to use to select the intrinsic dimensions.
#' @param init A character string or a vector of clusters. It is the way to initialize the E-M algorithm. There are five ways of initialization: \dQuote{kmeans} (default), \dQuote{param}, \dQuote{random}, \dQuote{mini-em} or \dQuote{vector}. See details for more information. It can also be directly initialized with a vector containing the prior classes of the observations.
#' @param init.vector A vector of integers or factors. It is a user-given initialization. It should be of the same length as of the data. Only used when \code{init="vector"}.
#' @param show Use show = FALSE to settle off the informations that may be printed.
#' @param mini.nb A vector of integers of length two. This parameter is used in the \dQuote{mini-em} initialization. The first integer sets how many times the algorithm is repeated; the second sets the maximum number of iterations the algorithm will do each time. For example, if init=\dQuote{mini-em} and mini.nb=c(5,10), the algorithm wil be lauched 5 times, doing each time 10 iterations; finally the algorithm will begin with the initialization that maximizes the log-likelihood. 
#' @param scaling Logical: whether to scale the dataset (mean=0 and standard-error=1 for each variable) or not. By default the data is not scaled.
#' @param min.individuals This parameter is used to control for the minimum population of a class. If the population of a class becomes stricly inferior to 'min.individuals' then the algorithm stops and gives the message: 'pop<min.indiv.'. Here the meaning of "population of a class" is the sum of its posterior probabilities. The value of 'min.individuals' cannot be lower than 2.
#' @param noise.ctrl This parameter avoids to have a too low value of the 'noise' parameter b. It garantees that the dimension selection process do not select too many dimensions (which leads to a potential too low value of the noise parameter b). When selecting the intrinsic dimensions using Cattell's scree-test or BIC, the function doesn't use the eigenvalues inferior to noise.ctrl, so that the intrinsic dimensions selected can't be higher or equal to the order of these eigenvalues.
#' @param mc.cores Positive integer, default is 1. If \code{mc.cores>1}, then parallel computing is used, using \code{mc.cores} cores. Warning for Windows users only: the parallel computing can sometimes be slower than using one single core (due to how parLapply works).
#' @param nb.rep A positive integer (default is 1). Each estimation (i.e. combination of (model, K, threshold)) is repeated \code{nb.rep} times and only the estimation with the highest log-likelihood is kept.
#' @param keepAllRes Logical. Should the results of all runs be kept? If so, an argument \code{all_results} is created in the results. Default is \code{TRUE}.
#' @param kmeans.control A list. The elements of this list should match the parameters of the kmeans initialization (see \code{\link[stats]{kmeans}} help for details). The parameters are \dQuote{iter.max}, \dQuote{nstart} and \dQuote{algorithm}.
#' @param d_max A positive integer. The maximum number of dimensions to be computed. Default is 100. It means that the instrinsic dimension of any cluster cannot be larger than \code{dmax}. It quickens a lot the algorithm for datasets with a large number of variables (e.g. thousands).
#' @param d DEPRECATED. This parameter is kept for retro compatibility. Now please use the parameter d_select. 
#' 
#' @details
#' Some information on the signification of the model names:
#' \describe{
#' 	\item{Akj are the parameters of the classes subspaces:}{
#' 		\itemize{
#' 		\item{if Akj: each class has its parameters and there is one parameter for each dimension}
#' 		\item{if Ak: the classes have different parameters but there is only one per class}
#' 		\item{if Aj: all the classes have the same parameters for each dimension (it's a particular case with a common orientation matrix)}
#' 		\item{if A: all classes have the same one parameter}
#' 		}
#' 	}
#' 
#' 	\item{Bk are the noises of the classes subspaces:}{
#' 		\itemize{
#' 			\item{If Bk: each class has its proper noise}
#' 			\item{if B:  all classes have the same noise}
#' 		}
#' 	}
#' 
#' 	\item{Qk is the orientation matrix of each class:}{ 
#' 		\itemize{
#' 			\item{if Qk: all classes have its proper orientation matrix}
#' 			\item{if Q: all classes have the same orientation matrix}
#' 		}
#' 	}
#' 	
#' 	\item{Dk is the intrinsic dimension of each class:}{ 
#' 		\itemize{
#' 			\item{if Dk: the dimensions are free and proper to each class}
#' 			\item{if D: the dimension is common to all classes}
#' 		}
#' 	}
#' }
#' The model \dQuote{ALL} will compute all the models, give their BIC and keep the model with the highest BIC value.
#' Instead of writing the model names, they can also be specified using an integer.  1 represents the most general model (\dQuote{AkjBkQkDk}) while 14 is the most constrained (\dQuote{ABQD}), the others  number/name matching are given below. Note also that several models can be run at once, by using a vector of models (e.g. model = c("AKBKQKD","AKJBQKDK","AJBQD") is equivalent to model = c(8,4,13); to run the 6 first models, use model=1:6). If all the models are to be run, model="all" is faster than model=1:14. 
#' \tabular{lcclc}{
#' AkjBkQkDk \tab   1   \tab   \tab  AkjBkQkD \tab   7   \cr 
#' AkBkQkDk \tab   2   \tab \tab  AkBkQkD \tab   8   \cr   
#' ABkQkDk \tab   3   \tab  \tab ABkQkD \tab   9   \cr   
#' AkjBQkDk \tab   4   \tab  \tab  AkjBQkD \tab   10   \cr   
#' AkBQkDk \tab   5   \tab  \tab  AkBQkD \tab   11   \cr   
#' ABQkDk \tab   6   \tab  \tab  ABQkD \tab   12  \cr
#' AjBQD \tab 13 \tab  \tab ABQD \tab 14
#' }
#' 
#' The parameter d, is used to select the intrinsic dimensions of the subclasses. Here are his definictions:
#' 		\itemize{
#' 			\item{\dQuote{Cattell}:}{
#' 				The Cattell's scree-test is used to gather the intrinsic dimension of each class. If the model is of common dimension (models 7 to 14), the scree-test is done on the covariance matrix of the whole dataset.
#' 			}
#' 			\item{\dQuote{BIC}:}{
#' 				The intrinsic dimensions are selected with the BIC criterion. See Bouveyron \emph{et al.} (2010) for a discussion of this topic.
#' 				For common dimension models, the procedure is done on the covariance matrix of the whole dataset.
#' 			}
#' 			\item{Note that "Cattell" (resp. "BIC") can be abreviated to "C" (resp. "B") and that this argument is not case sensitive.}
#' 		}
#' 		
#' The different initializations are:
#' \describe{
#' 	\item{\dQuote{param}:}{it is initialized with the parameters, the means being generated by a multivariate normal distribution and the covariance matrix being common to the whole sample}
#' 
#' 	\item{\dQuote{mini-em}:}{it is an initialization strategy, the classes are randomly initialized and the E-M algorithm makes several iterations, this action is repetead a few times (the default is 5 iterations and 10 times), at the end, the initialization choosen is the one which maximise the log-likelihood (see mini.nb for more information about its parametrization)}
#' 
#' 	\item{\dQuote{random}:}{the classes are randomly given using a multinomial distribution}
#' 
#' 	\item{\dQuote{kmeans}:}{the classes are initialized using the kmeans function (with: algorithm="Hartigan-Wong"; nstart=4; iter.max=50); note that the user can use his own arguments for kmeans using the dot-dot-dot argument } 
#' 	
#' 	\item{A prior class vector:}{It can also be directly initialized with a vector containing the prior classes of the observations. To do so use \code{init="vector"} and provide the vector in the argument \code{init.vector}.}
#' }
#' 
#' The BIC criterion used in this function is to be maximized and is defined as 2*LL-k*log(n) where LL is the log-likelihood, k is the number of parameters and n is the number of observations.
#' 
#'
#' @references
#' Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519
#' 
#' Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, \bold{46}(6), 1--29, url: \href{http://www.jstatsoft.org/v46/i06/}{http://www.jstatsoft.org/v46/i06/}
#' 
#' @author
#' Laurent Berge, Charles Bouveyron and Stephane Girard
#' 
#' @seealso
#' \code{\link{hdda}}, \code{\link{predict.hdc}}, \code{\link{plot.hdc}}.
#' 
#' @return
#' hddc returns an 'hdc' object; it's a list containing:
#' \item{ model }{The name of the model.}
#' \item{ K }{The number of classes.}
#' \item{ d }{The dimensions of each class.}
#' \item{ a }{The parameters of each class subspace.}
#' \item{ b }{The noise of each class subspace.}
#' \item{ mu }{The mean of each variable for each class.}
#' \item{ prop }{The proportion of each class.}
#' \item{ ev }{The eigen values of the var/covar matrix.}
#' \item{ Q }{The orthogonal matrix of orientation of each class.}
#' \item{ loglik }{The log-likelihood.}
#' \item{ loglik_all }{The log-likelihood of all iterations.}
#' \item{ posterior }{The matrix of the probabilities to belong to a class for each observation and each class.}
#' \item{ class }{The class vector obtained by the clustering.}
#' \item{ com_ev }{Only if this is a common dimension model. The eigenvalues of the var/covar matrix of the whole dataset.}
#' \item{ N }{The number of observations.}
#' \item{ complexity }{The number of parameters of the model.}
#' \item{ threshold }{The threshold used for the Cattell scree-test.}
#' \item{ d_select }{The way the dimensions were selected.}
#' \item{ BIC }{The BIC of the model.}
#' \item{ ICL }{The ICL of the model.}
#' \item{ criterion }{The criterion used to select the model.}
#' \item{ call }{The call.}
#' \item{ complexity_allModels }{The number of parameters for all the models that were run (used to compute the slope heuristic).}
#' \item{ allCriteria }{The data.frame with the combination (model, K, threshold) and the associated values of the likelihood (LL), BIC and ICL, as well as the rank of each of the models with respect to the selection criterion.}
#' \item{ all_results }{Only if \code{keepAllRes=TRUE}. The parameters of all estimations that were run.}
#' \item{ scaling }{Only if \code{scaling=TRUE}. The centers and the standard deviation of the original dataset.}
#' 
#' @concept
#' clustering
#'
#' @examples
#' # Example 1:
#' data <- simuldata(1000, 1000, 50)
#' X <- data$X
#' clx <- data$clx
#' Y <- data$Y
#' cly <- data$cly
#' 
#' #clustering of the simulated dataset:
#' prms1 <- hddc(X, K=3, algo="CEM", init='param')                
#' 
#' #class vector obtained by the clustering:
#' prms1$class                
#' 
#' #We can look at the adjusted rand index to assess the goodness of fit
#' res1 <- predict(prms1, X, clx)
#' res2 <- predict(prms1, Y)       
#' #the class predicted using hddc parameters on the test dataset:  
#' res2$class                                                           
#' 
#' 
#' # Example 2:
#' data(Crabs)
#' 
#' # clustering of the Crabs dataset:
#' prms3 <- hddc(Crabs[,-1], K=4, algo="EM", init='mini-em')        
#' res3 <- predict(prms3, Crabs[,-1], Crabs[,1])
#' 
#' # another example using the Crabs dataset
#' prms4 <- hddc(Crabs[,-1], K=1:8, model=c(1,2,7,9))
#' 
#' # model=c(1,2,7,9) is equivalent to:
#' # model=c("AKJBKQKDK","AKBKQKDK","AKJBKQKD"#' ,"ABKQKD") 
#' res4 <- predict(prms4, Crabs[,-1], Crabs[,1])
#' 
#' # PARALLEL COMPUTING
#' \dontrun{
#' # Same example but with Parallel Computing => platform specific
#' # (slower for Windows users)
#' # To enable it, just use the argument 'mc.cores'
#' prms5 <- hddc(Crabs[,-1], K=1:8, model=c(1,2,7,9), mc.cores=2)
#' }
#' 
hddc  <- function(data, K=1:10, model=c("AkjBkQkDk"), threshold=0.2, criterion="bic", com_dim=NULL, itermax=200, eps=1e-3, algo='EM', d_select="Cattell", init='kmeans', init.vector, show=TRUE, mini.nb=c(5, 10), scaling=FALSE, min.individuals=2, noise.ctrl=1e-8, mc.cores=1, nb.rep=1, keepAllRes=TRUE, kmeans.control = list(), d_max=100, d){
	
	# For compatibility with old versions of HDclassif
	if(!missing(d) & missing(d_select)) d_select = d
	
	#
	# CONTROLS
	#
	
	call = match.call()
	hddc_control(call)
	# Control of match.args:
	criterion = myAlerts(criterion, "criterion", "singleCharacterMatch.arg", "HDDC: ", c("bic", "icl"))
	algo = myAlerts(algo, "algo", "singleCharacterMatch.arg", "HDDC: ", c('EM', 'CEM', 'SEM'))
	d_select = myAlerts(d_select, "d_select", "singleCharacterMatch.arg", "HDDC: ", c("cattell", "bic"))
	init = myAlerts(init, "init", "singleCharacterMatch.arg", "HDDC: ", c('random', 'kmeans', 'mini-em', 'param', "vector"))
	# We get the model names, properly ordered
	model = hdc_getTheModel(model, all2models = TRUE)
	# kmeans controls
	kmeans.control = default_kmeans_control(kmeans.control)
	
	data <- as.matrix(data)
	if (scaling) {
		data <- scale(data)
		scaling <- list(mu=attr(data, "scaled:center"), sd=attr(data, "scaled:scale"))
	} else scaling <- NULL
	
	BIC <- ICL <- c()
	p <- ncol(data)
	
	#
	# Preparing the parallel
	#
	
	if(d_select=="bic"){
		# If the dimension selection is done with BIC, we don't care of the threshold
		threshold = "bic"
	} 
	
	if(max(table(K))>1) warning("The number of clusters, K, is made unique (repeated values are not tolerated).")
	K = sort(unique(K))
	if(any(K==1)){
		# K=1 => only one model required
		K = K[K!=1]
		addrows = data.frame(model="AKJBKQKDK", K=1, threshold)
	} else {
		addrows = c()
	}
	
	mkt_Expand = expand.grid(model=model, K=K, threshold=threshold)
	mkt_Expand = do.call(rbind, replicate(nb.rep, mkt_Expand, simplify=FALSE))
	mkt_Expand = rbind(addrows, mkt_Expand) #no need for several runs for K==1
	
	model = as.character(mkt_Expand$model)
	K = mkt_Expand$K
	threshold = mkt_Expand$threshold
	
	# We transform it into an univariate form
	mkt_univariate = apply(mkt_Expand, 1, paste, collapse= "_")
	
	# Mon 'caller' for LTBM that will be used in multi-cores lapply
	hddcWrapper = function(mkt_univariate, ...){
		
		mkt_splitted =  strsplit(mkt_univariate, "_")
		
		# on retrouve model, K and threshold
		model = sapply(mkt_splitted, function(x) x[1])
		K = sapply(mkt_splitted, function(x) as.numeric(x[2]))
		threshold = sapply(mkt_splitted, function(x) ifelse(x[3]=="bic","bic",as.numeric(x[3])))
		
		# (::: is needed for windows multicore)
		res = "unknown error"
		# try(res <- HDclassif:::hddc_main(model=model, K=K, threshold=threshold, ...))
		try(res <- hddc_main(model=model, K=K, threshold=threshold, ...))
		res
	}
	
	# We reset the number of cores to use
	nRuns = length(mkt_univariate)
	if(nRuns<mc.cores) mc.cores = nRuns
	
	# We swicth to the right number of cores + a warning if necessary
	max_nb_of_cores = parallel::detectCores()
	if(mc.cores>max_nb_of_cores){
		warning("The argument mc.cores is greater than its maximun.\nmc.cores was set to ", max_nb_of_cores)
		mc.cores = max_nb_of_cores
	}
	
	
	#
	# Parallel estimations
	#
	
	if(mc.cores == 1){
		# If there is no need for parallel, we just use lapply // in order to have the same output
		
		par.output = lapply(mkt_univariate, hddcWrapper, DATA=data, method=d_select, algo=algo, itermax=itermax, eps=eps, init=init, init.vector=init.vector, mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max)
		
	} else if(Sys.info()[['sysname']] == 'Windows'){
		# we use parLapply:
		
		## create clusters
		cl = parallel::makeCluster(mc.cores)
		## load the packages
		loadMyPackages = function(x){
			# loadMyPackages = function(none, myFuns){
			# we load the package
			library(HDclassif)
			# we add the functions in the global env
			# for(i in 1:length(myFuns)){
			# 	funName = names(myFuns)[i]
			# 	fun = myFuns[[i]]
			# 	assign(funName, fun, .GlobalEnv)
			# }
		}
		## Create the functions to export
		# myFuns = list(hddc_main = hddc_main, hddc_e_step=hddc_e_step, hddc_m_step=hddc_m_step, hdc_getComplexity=hdc_getComplexity, hdc_myEigen=hdc_myEigen, hdclassif_dim_choice=hdclassif_dim_choice, hdclassif_bic=hdclassif_bic)
		# par.setup = parallel::parLapply(cl, 1:length(cl), loadMyPackages, myFuns=myFuns)
		par.setup = parallel::parLapply(cl, 1:length(cl), loadMyPackages)
		
		## run the parallel
		par.output = NULL
		try(par.output <- parallel::parLapply(cl, mkt_univariate, hddcWrapper, DATA=data, method=d_select, algo=algo, itermax=itermax, eps=eps, init=init, init.vector=ifelse(missing(init.vector), NA, init.vector), mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max))
		
		## Stop the clusters
		parallel::stopCluster(cl)
		
		if(is.null(par.output)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
		
	} else {
		# we use mclapply
		
		par.output = NULL
		try(par.output <- parallel::mclapply(mkt_univariate, hddcWrapper, DATA=data, method=d_select, algo=algo, itermax=itermax, eps=eps, init=init, init.vector=init.vector, mini.nb=mini.nb, min.individuals=min.individuals, noise.ctrl=noise.ctrl, com_dim=com_dim, kmeans.control=kmeans.control, d_max=d_max, mc.cores=mc.cores))
		
		if(is.null(par.output)) stop("Unknown error in the parallel computing. Try mc.cores=1 to detect the problem.")
		
	}
	
	#
	# The results are retrieved
	#
	
	getElement = function(x, what, valueIfNull = -Inf){
		# attention si x est le modele nul
		if(length(x)==1) return(valueIfNull)
		if(!is.list(x) && !what %in% names(x)) return(NA)
		x[[what]][length(x[[what]])]
	}
	
	getComment = function(x){
		# we get the error message
		if(length(x)==1) return(x)
		return("")
	}
	
	# All likelihoods
	LL_all = sapply(par.output, getElement, what="loglik")
	comment_all = sapply(par.output, getComment)
	
	# If no model is valid => problem
	if(all(!is.finite(LL_all))){
		warning("All models diverged.")
		allCriteria = data.frame(model=model, K=K, threshold=threshold, LL = LL_all, BIC=NA, comment=comment_all)
		res = list()
		res$allCriteria = allCriteria
		return(res)
	}
	
	# We select, for each (Q,K), the best run
	n = nrow(mkt_Expand)
	modelKeep = sapply(unique(mkt_univariate), function(x) (1:n)[mkt_univariate==x][which.max(LL_all[mkt_univariate==x])])
	# => we select only the best models
	LL_all = LL_all[modelKeep]
	comment_all = comment_all[modelKeep]
	par.output = par.output[modelKeep]
	BIC = sapply(par.output, getElement, what="BIC")
	ICL = sapply(par.output, getElement, what="ICL")
	comp_all = sapply(par.output, getElement, what="complexity", valueIfNull=NA)
	model = model[modelKeep]
	threshold = threshold[modelKeep]
	K = K[modelKeep]
	
	# We define the criterion of model selection
	CRIT = switch(criterion,
				  bic = BIC,
				  icl = ICL)
	
	# The order of the results
	myOrder = order(CRIT, decreasing = TRUE)
	
	# On sauvegarde le bon modele + creation de l'output
	qui = which.max(CRIT)
	
	prms = par.output[[qui]]
	prms$criterion = CRIT[qui]
	names(prms$criterion) = criterion
	
	# Other output
	prms$call = call
	# We add the complexity
	names(comp_all) = mkt_univariate[modelKeep]
	prms$complexity_allModels = comp_all
	
	# Display
	if(show){
		if(n>1) cat("HDDC: \n")
		
		model2print = sapply(model, function(x) sprintf("%*s", max(nchar(model)), x))
		K2print = as.character(K)
		K2print = sapply(K2print, function(x) sprintf("%*s", max(nchar(K2print)), x))
		thresh2print = as.character(threshold)
		thresh_width = max(nchar(thresh2print))
		thresh2print = sapply(thresh2print, function(x) sprintf("%s%s", x, paste0(rep("0", thresh_width - nchar(x)), collapse="") ))
		
		# on cree une data.frame
		myResMat = cbind(model2print[myOrder], K2print[myOrder], thresh2print[myOrder], addCommas(CRIT[myOrder]), comment_all[myOrder])
		
		myResMat = as.data.frame(myResMat)
		names(myResMat) = c("model", "K", "threshold", toupper(criterion), "comment")
		row.names(myResMat) = 1:nrow(myResMat)
		
		# if no problem => no comment
		if(all(comment_all == "")) myResMat$comment = NULL
		
		print(myResMat)
		
		msg = switch(criterion, bic="BIC", icl="ICL")
		cat("\nSELECTED: model ", prms$model, " with ", prms$K, " clusters.\n")
		cat("Selection Criterion: ", msg, ".\n", sep="")
		
	}
	
	# We also add the matrix of all criteria
	allCriteria = data.frame(model=model[myOrder], K=K[myOrder], threshold=threshold[myOrder], LL=LL_all[myOrder], BIC=BIC[myOrder], ICL=ICL[myOrder], rank = 1:length(myOrder))
	
	# we add the comments if necessary
	if(any(comment_all != "")) allCriteria$comment = comment_all[myOrder]
	prms$allCriteria = allCriteria
	
	# If all the results are kept
	if(keepAllRes){
		all_results = par.output
		names(all_results) = mkt_univariate[modelKeep]
		prms$all_results = all_results
	}
	
	# Other stuff
	prms$scaling <- scaling
	prms$threshold <- threshold[qui]
	
	return(prms)
}

hddc_main <- function(DATA, K, model, threshold, method, algo, itermax, eps, init, init.vector, mini.nb, min.individuals, noise.ctrl, com_dim=NULL, kmeans.control, d_max, ...){ 
	
	# for debug
	debug = FALSE

	ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD")
	p <- ncol(DATA)
	N <- nrow(DATA)
	com_ev <- NULL
	
	# We set d_max to a proper value
	d_max = min(N, p, d_max)
	
	if ( any(model==ModelNames[7:14]) ){
		# Common dimension models
		MU <- colMeans(DATA)
		if (N<p) {
			Y <- (DATA-matrix(MU, N, p, byrow=TRUE))/sqrt(N)
			YYt <- tcrossprod(Y)
			com_ev <- hdc_myEigen(YYt, d_max, only.values = TRUE)$values
		} else{
			S <- crossprod(DATA-matrix(MU, N, p, byrow=TRUE))/N
			com_ev <- hdc_myEigen(S, d_max, only.values = TRUE)$values
		}
		if(is.null(com_dim)) com_dim <- hdclassif_dim_choice(com_ev, N, method, threshold, FALSE, noise.ctrl)
	}
	
	if (K>1){
		t <- matrix(0, N, K)
		if(init == "vector"){
			init.vector = unclass(init.vector)
			name <- unique(init.vector)
			for (i in 1:K) t[which(init.vector==name[i]), i] <- 1
		} else if (init=='param'){
			MU <- colMeans(DATA)
			prop <- rep(1/K, K)
			S <- crossprod(DATA - matrix(MU, N, p, byrow=TRUE))/N
			donnees <- eigen(S, symmetric=TRUE)
			ev <- donnees$values
			d <- if(is.numeric(method)) method else hdclassif_dim_choice(ev, N, method, threshold, FALSE, noise.ctrl)
			a <- ev[1:d]
			b <- sum(ev[(d[1]+1):p])/(p-d[1])
			
			Q <- donnees$vectors[, 1:d]
			mu <- mvrnorm(K, MU, S)
			
			K_pen <- diag((mu%*%Q%*%diag(1/a, d, d))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a, d, d))%*%(t(Q)%*%t(DATA))+1/b*(diag(tcrossprod(mu))-2*mu%*%t(DATA)+2*(mu%*%Q)%*%(t(Q)%*%t(DATA))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
			
			t <- matrix(0, N, K)
			for (i in 1:K) t[, i]=1/rowSums(exp((K_pen[i, ]-t(K_pen))/2))
		} else if (init=='kmeans') {
			kmc = kmeans.control
			cluster <- kmeans(DATA, K, iter.max=kmc$iter.max, nstart=kmc$nstart, algorithm=kmc$algorithm, trace=kmc$trace)$cluster
			for (i in 1:K) t[which(cluster==i), i] <- 1
		} else if (init=='mini-em'){
			prms_best <- 1
			for (i in 1:mini.nb[1]){
				prms <- hddc_main(DATA, K, model, threshold, method, algo, mini.nb[2], 0, 'random', mini.nb = mini.nb, min.individuals = min.individuals, noise.ctrl = noise.ctrl, com_dim = com_dim, d_max=d_max)
				if(length(prms)!=1){
					if (length(prms_best)==1) prms_best <- prms
					else if (prms_best$loglik[length(prms_best$loglik)]<prms$loglik[length(prms$loglik)]) prms_best <- prms
				}
			}
			
			if (length(prms_best)==1) return(1)
			t <- prms_best$posterior
		} else {
			t <- t(rmultinom(N, 1, rep(1/K, K)))
			compteur=1
			while(min(colSums(t))<1 && (compteur <- compteur+1)<5) t <- t(rmultinom(N, 1, rep(1/K, K)))
			if(min(colSums(t))<1) return("Random initialization failed (n too small)")
		}
	} else t <- matrix(1, N, 1)
	
	likely <- c()
	I <- 0
	test <- Inf
	while ((I <- I+1)<=itermax && test>=eps){
		
		if (algo!='EM' && I!=1) t <- t2
		
		if(debug) cat("Cluster sizes: ", colSums(t), "\n")
		
		# Error catching
		if (K>1){
			if(any(is.na(t))) return("unknown error: NA in t_ik")
			
			if(any(colSums(t>1/K)<min.individuals)) return("pop<min.individuals")
		}
		
		m <- hddc_m_step(DATA, K, t, model, threshold, method, noise.ctrl, com_dim, d_max)
		t <- hddc_e_step(DATA, m)
		L <- t$L
		t <- t$t
		if (algo=='CEM') {
			t2 <- matrix(0, N, K)
			t2[cbind(1:N, max.col(t))] <- 1
		} else if(algo=='SEM') { 
			t2 <- matrix(0, N, K)
			for (i in 1:N)	t2[i, ] <- t(rmultinom(1, 1, t[i, ]))
		}
		likely[I] <- L
		if (I!=1) test <- abs(likely[I]-likely[I-1])
		
		if(debug){
			print("d=")
			print(m$d)
			print("a=")
			print(m$a)
			print("b=")
			print(m$b)
		}
		
	}
	
	# We retrieve the parameters
	
	# a
	if ( model%in%c('AKBKQKDK', 'AKBQKDK', 'AKBKQKD', 'AKBQKD') ) {
		a <- matrix(m$a[, 1], 1, m$K, dimnames=list(c("Ak:"), 1:m$K))
	} else if(model=='AJBQD') {
		a <- matrix(m$a[1, ], 1, m$d[1], dimnames=list(c('Aj:'), paste('a', 1:m$d[1], sep='')))
	} else if ( model%in%c('ABKQKDK', 'ABQKDK', 'ABKQKD', 'ABQKD', "ABQD") ) {
		a <- matrix(m$a[1], dimnames=list(c('A:'), c('')))
	} else a <- matrix(m$a, m$K, max(m$d), dimnames=list('Class'=1:m$K, paste('a', 1:max(m$d), sep='')))
	
	# b
	if ( model%in%c('AKJBQKDK', 'AKBQKDK', 'ABQKDK', 'AKJBQKD', 'AKBQKD', 'ABQKD', 'AJBQD', "ABQD") ) {
		b <- matrix(m$b[1], dimnames=list(c('B:'), c('')))
	} else b <- matrix(m$b, 1, m$K, dimnames=list(c("Bk:"), 1:m$K))
	
	# d, mu, prop
	d <- matrix(m$d, 1, m$K, dimnames=list(c('dim:'), "Intrinsic dimensions of the classes:"=1:m$K))
	mu <- matrix(m$mu, m$K, p, dimnames=list('Class'=1:m$K, 'Posterior group means:'=paste('V', 1:p, sep='')))
	prop <- matrix(m$prop, 1, m$K, dimnames=list(c(''), 'Posterior probabilities of groups'=1:m$K))
	
	# Other elements
	complexity <- hdc_getComplexity(m, p)
	class(b) <- class(a) <- class(d) <- class(prop) <- class(mu) <- 'hd'
	cls <- max.col(t)
	converged = test<eps
	
	params = list(model=model, K=K, d=d, a=a, b=b, mu=mu, prop=prop, ev=m$ev, Q=m$Q, loglik=likely[length(likely)], loglik_all = likely, posterior=t, class=cls, com_ev=com_ev, N=N, complexity=complexity, threshold=threshold, d_select=method, converged=converged)

	# We compute the BIC / ICL
	bic_icl = hdclassif_bic(params, p, DATA)
	params$BIC = bic_icl$bic
	params$ICL = bic_icl$icl
	
	# We set the class
	class(params) <- 'hdc'
	
	return(params)
}

hddc_e_step  <- function(x, par){
	p <- ncol(x)
	N <- nrow(x)
	K <- par$K
	a <- par$a
	b <- par$b
	mu <- par$mu
	d <- par$d
	prop <- par$prop
	Q <- par$Q
	
	b[b<1e-6] <- 1e-6
	
	if(par$model=="AJBQD") {
		
		K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(x))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)+2*(mu%*%Q)%*%(t(Q)%*%t(x))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		
	} else if(par$model=="ABQD") {
		
		K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(x)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(x))
		
	} else{
		K_pen <- matrix(0,K,N)
		for (i in 1:K) {
			s <- sum(log(a[i,1:d[i]]))
			X <- x - matrix(mu[i,], N, p, byrow=TRUE)
			proj <- (X%*%Q[[i]])%*%t(Q[[i]])
			A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
			B <- X-proj
			K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
		}
	}
	
	A <- -1/2*t(K_pen)
	L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	
	t <- matrix(0,N,K)
	for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
	list(t=t,L=L)
}

hddc_m_step  <- function(x, K, t, model, threshold, method, noise.ctrl, com_dim, d_max){

	# Some parameters
	N <- nrow(x)
	p <- ncol(x)
	prop <- c()
	n <- colSums(t)
	prop <- n/N
	mu <- matrix(NA, K, p)
	for (i in 1:K) mu[i, ] <- colSums(x*t[, i])/n[i]
	
	ind <- apply(t>0, 2, which)
	n_bis <- c()
	for(i in 1:K) n_bis[i] <- length(ind[[i]])
	
	#
	# Calculation on Var/Covar matrices
	#
	
	# we keep track of the trace (== sum of eigenvalues) to compute the b
	traceVect = c()
	
	if (N<p) {
		if( model%in%c("AJBQD", "ABQD") ){
			Y <- matrix(0, N, p)
			for (i in 1:K) Y <- Y+(x-matrix(mu[i, ], N, p, byrow=TRUE))/sqrt(N)*sqrt(t[, i])
			YYt = tcrossprod(Y)
			donnees <- hdc_myEigen(YYt, d_max)
			traceVect = sum(diag(YYt))
			ev <- donnees$values
		} else {
			Y <- vector(mode='list', length=K)
			# ev <- matrix(0, K, N) # now we use d_max
			ev <- matrix(0, K, d_max)
			Q <- vector(mode='list', length=K)
			for (i in 1:K){ 
				Y[[i]] <- (x-matrix(mu[i, ], N, p, byrow=TRUE))/sqrt(n[i])*sqrt(t[, i])
				YYt = tcrossprod(Y[[i]])
				donnees <- hdc_myEigen(YYt, d_max)
				traceVect[i] = sum(diag(YYt))
				# ev[i, 1:N] <- donnees$values # now we use d_max
				ev[i, ] <- donnees$values
				Q[[i]] <- donnees$vectors
			}
		}
	} else if ( model%in%c("AJBQD", "ABQD") ){
		W <- matrix(0, p, p)
		for (i in 1:K) W <- W + crossprod((x-matrix(mu[i, ], N, p, byrow=TRUE))*sqrt(t[, i]))/N
		donnees <- hdc_myEigen(W, d_max)
		traceVect = sum(diag(W))
		ev <- donnees$values
	} else {
		# ev <- matrix(0, K, p) # now we use d_max
		ev <- matrix(0, K, d_max)
		Q <- vector(mode='list', length=K)
		for (i in 1:K){ 
			W = crossprod((x-matrix(mu[i, ], N, p, byrow=TRUE))*sqrt(t[, i]))/n[i]
			# browser()
			donnees <- hdc_myEigen(W, d_max)
			traceVect[i] = sum(diag(W))
			ev[i, ] <- donnees$values
			Q[[i]] <- donnees$vectors
		}
	}	
	
	# Intrinsic dimensions selection
	# browser()
	if (model%in%c("AJBQD", "ABQD")){
		d <- rep(com_dim, length=K)
	} else if ( model%in%c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", "AKBQKD", "ABQKD") ){
		dmax <- min(apply((ev>noise.ctrl)*rep(1:ncol(ev), each=K), 1, which.max))-1
		if(com_dim>dmax) com_dim <- max(dmax, 1)
		d <- rep(com_dim, length=K)
	} else {
		d <- hdclassif_dim_choice(ev, n, method, threshold, FALSE, noise.ctrl)
	}
	
	# Setup of the Qi matrices	
	
	if ( model%in%c("AJBQD", "ABQD") ){
		if (N>=p) Q <- matrix(donnees$vectors[, 1:d[1]], p, d[1])
		else {
			Q <- matrix(t(Y)%*%donnees$vectors[, 1:d[1]], p, d[1])
			normalise <- c()
			for(i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[, i]))
			Q <- Q/matrix(sqrt(normalise), p, d, byrow=TRUE)
		}
	} else if(N>=p) {
		for(i in 1:K) Q[[i]] <- matrix(Q[[i]][, 1:d[i]], p, d[i])
	} else{
		for (i in 1:K){ 
			Q[[i]] <- t(Y[[i]])%*%(Q[[i]][, 1:d[i]])
			normalise <- c()
			for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][, j])))
			Q[[i]] <- Q[[i]]/matrix(sqrt(normalise), p, d[i], byrow=TRUE)
		}
	}
	
	# Calculation of the remaining parameters of the selected model	
	
	# PARAMETER a
	ai <- matrix(NA, K, max(d))
	if ( model%in%c('AKJBKQKDK', 'AKJBQKDK', 'AKJBKQKD', 'AKJBQKD') ){
		for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
	} else if ( model%in%c('AKBKQKDK', 'AKBQKDK' , 'AKBKQKD', 'AKBQKD') ){
		for (i in 1:K) ai[i, ] <- rep(sum(ev[i, 1:d[i]])/d[i], length=max(d))
	} else if(model=="AJBQD"){
		for (i in 1:K) ai[i, ] <- ev[1:d[1]]
	} else if(model=="ABQD") {
		ai[] <- sum(ev[1:d[1]])/d[1]
	} else {
		a <- 0
		eps <- sum(prop*d)
		for (i in 1:K) a <- a + sum(ev[i, 1:d[i]])*prop[i]
		ai <- matrix(a/eps, K, max(d))
	}
	
	# PARAMETER b
	bi <- c()
	# denom = min(N,p) # DEPREC => ask Charles
	denom = p
	if ( model%in%c('AKJBKQKDK', 'AKBKQKDK', 'ABKQKDK', 'AKJBKQKD', 'AKBKQKD', 'ABKQKD') ){
		for(i in 1:K){
			remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
			# bi[i] <- sum(ev[i, (d[i]+1):min(N, p)])/(p-d[i])
			bi[i] <- remainEV/(denom-d[i])
		} 
	} else if ( model%in%c("ABQD", "AJBQD") ){
		remainEV = traceVect - sum(ev[1:d[1]])
		# bi[1:K] <- sum(ev[(d[1]+1):min(N, p)])/(min(N, p)-d[1])
		bi[1:K] <- remainEV/(denom-d[1])
	} else {		
		b <- 0
		eps <- sum(prop*d)
		for(i in 1:K){
			remainEV = traceVect[i] - sum(ev[i, 1:d[i]])
			# b <- b + sum(ev[i, (d[i]+1):min(N, p)])*prop[i]
			b <- b + remainEV*prop[i]
		}
		bi[1:K] <- b/(denom-eps)
	}
	
	# We adjust the values of b if they are too low
	bi[bi<noise.ctrl] = noise.ctrl
	
	list(model=model, K=K, d=d, a=ai, b=bi, mu=mu, prop=prop, ev=ev, Q=Q)
}

hddc_ari <- function(x,y){
	#This function is drawn from the mclust package
	x <- as.vector(x)
	y <- as.vector(y)
	tab <- table(x, y)
	if (all(dim(tab) == c(1, 1))) return(1)
	a <- sum(choose(tab, 2))
	b <- sum(choose(rowSums(tab), 2)) - a
	c <- sum(choose(colSums(tab), 2)) - a
	d <- choose(sum(tab), 2) - a - b - c
	ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
	return(ARI)
}

#' Slope Heuristic for HDDC objects
#' 
#' This function computes the slope heuristic for an 
#'
#' @param x An \code{hdc} object, obtained from the function \code{\link{hddc}}.
#'
#' @return
#' A list of two elements:
#' \item{best_model}{The index of the best model, among all estimated models.}
#' \item{crit}{The slope heuristic criterion for each of the estimated models.}
#'  
#' @export
#'
#' @examples
#' # Clustering of the Crabs data set
#' data(Crabs)
#' prms = hddc(Crabs[,-1], K=1:10)
#' slope = slopeHeuristic(prms)
#' plot(slope$crit) # The best model is indeed for 4 clusters
#' prms$all_results[[slope$best_model]] # we extract the best model
#' 
#' 
slopeHeuristic <- function(x) {
	
	# x is a hddc object
	
	nbparam = x$complexity_allModels
	loglik = x$allCriteria$LL
	
	# Controlling for NA
	n = length(nbparam)
	crit = rep(NA, n)
	quiNotNA = which(!is.na(nbparam))
	
	nbparam = nbparam[quiNotNA]
	loglik = loglik[quiNotNA]
	
	# Slope heuristic
	dd = data.frame(nbp=nbparam,ll=loglik)
	fit = MASS::rlm(ll ~ nbp, data=dd, method='MM')
	if(fit$coefficients[2]<0) fit$coefficients[2]=0
	dd$llpen = dd$ll- 2* fit$coefficients[2]*dd$nbp
	
	crit[quiNotNA] = 2*dd$llpen
	
	names(crit) = names(nbparam)
	
	list(best_model = which.max(crit), crit = crit)
}

####
#### CONTROLS ####
####

hddc_control = function(call){
	
	prefix = "HDDC: "
	myCallAlerts(call, "data", "matrix,data.frame", 3, TRUE, prefix)
	myCallAlerts(call, "K", "integerVector", 3, FALSE, prefix)
	myCallAlerts(call, "model", "vector", 3, FALSE, prefix)
	myCallAlerts(call, "threshold", "numericVectorGE0LE1", 3, FALSE, prefix)
	myCallAlerts(call, "criterion", "character", 3, FALSE, prefix)
	myCallAlerts(call, "com_dim", "singleIntegerGE1", 3, FALSE, prefix)
	myCallAlerts(call, "itermax", "singleIntegerGE0", 3, FALSE, prefix)
	myCallAlerts(call, "eps", "singleNumericGE0", 3, FALSE, prefix)
	myCallAlerts(call, "graph", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "algo", "singleCharacter", 3, FALSE, prefix)
	myCallAlerts(call, "d_select", "singleCharacter", 3, FALSE, prefix)
	myCallAlerts(call, "init", "singleCharacter", 3, FALSE, prefix)
	myCallAlerts(call, "show", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "mini.nb", "integerVectorGE1", 3, FALSE, prefix)
	myCallAlerts(call, "scaling", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "min.individuals", "singleIntegerGE2", 3, FALSE, prefix)
	myCallAlerts(call, "noise.ctrl", "singleNumericGE0", 3, FALSE, prefix)
	myCallAlerts(call, "mc.cores", "singleIntegerGE1", 3, FALSE, prefix)
	myCallAlerts(call, "nb.rep", "singleIntegerGE1", 3, FALSE, prefix)
	myCallAlerts(call, "keepAllRes", "singleLogical", 3, FALSE, prefix)
	myCallAlerts(call, "d_max", "singleIntegerGE1", 3, FALSE, prefix)
	
	
	####
	#### SPECIFIC controls 
	####
	
	# Getting some elements
	data = eval.parent(call[["data"]], 2)
	K = eval.parent(call[["K"]], 2)
	init = eval.parent(call[["init"]], 2)
	criterion = eval.parent(call[["criterion"]], 2)
	
	# No NA in the data:
	if (any(is.na(data))) stop("NA values in the data are not supported. Please remove them beforehand.")
	
	# Size of the data
	if(any(K>2*NROW(data))) stop("The number of observations must be at least twice the number of clusters ")

	# Initialization Controls
	if(!is.null(init)){
		
		# we get the value of the initialization
		init = myAlerts(init, "init", "singleCharacterMatch.arg", "HDDC: ", c('random', 'kmeans', 'mini-em', 'param', "vector"))
		
		# Custom initialization => controls and setup
		if(init == "vector"){
			myCallAlerts(call, "init.vector", "(integer,factor)Vector", 3, FALSE, prefix)
			
			init.vector = eval.parent(call[["init.vector"]], 2)
			
			if(is.null(init.vector)) stop("HDDC: When init='vector', the argument 'init.vector' should be provided.")
			
			if(length(unique(K))>1) stop("HDDC: Several number of classes K cannot be estimated when init='vector'.")
			
			init.vector <- unclass(as.factor(init.vector))
			if(K!=max(init.vector)) stop("The number of class K, and the number of classes in the initialization vector are different")
			
			if( length(init.vector)!=nrow(data) ) stop("The size of the initialization vector is different of the size of the data")
		}
		
		# The param init
		if (init=='param' && nrow(data)<ncol(data)){
			stop("The 'param' initialization can't be done when N<p")
		}
		
		# The mini.em init
		if (init=='mini-em'){
			
			mini.nb = eval.parent(call[["mini.nb"]], 2)
			
			if(!is.null(mini.nb) && length(mini.nb)!=2){
				stop("The parameter mini.nb must be a vector of length 2 with integers\n")
			}
			
		}
	}
	
}

default_kmeans_control = function(control){
	
	myAlerts(control,"kmeans.control","list","kmeans controls: ")
	
	#
	# Default values of the control parameters
	#
	
	myDefault = list()
	myDefault$iter.max = 10
	myDefault$nstart = 1
	myDefault$algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen")
	myDefault$trace = FALSE
	
	#
	# Types of each arg
	#
	
	myTypes = c("singleIntegerGE1", "singleIntegerGE1", "match.arg", "singleLogical")
	
	#
	# Recreation of the kmeans controls + Alerts
	#
	
	control = matchTypeAndSetDefault(control, myDefault, myTypes, "kmeans list of controls: ")
	
	return(control)
}
