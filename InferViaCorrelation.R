#dyn.load("InferViaCorrelation.so")

InferViaCorrelation <- function(data,treshold=0)
{	# data is one timeseries as created by "Boolnet:generateTimeSeries[[1]]"
	# requires named matrix
	if (is.matrix(data_input))
	{	data <- list(data_input)
	}
	else
	{	data <- data_input
	}
	if (is.null(rownames(data)))
	{	stop("rownames required")
	}
	result <- .Call("InferViaCorrelation",
		data,
		dim(data)[1],
		dim(data)[2],
		rownames(data),
		treshold
		)
	return(toBoolNet(names(result),result))
}

InferViaCorrelation_multiple_ts <- function(data,treshold=0,treshold_multi=0.5)
{	# data is timeseries as created by "Boolnet:generateTimeSeries"
	# requires list of named matrices
	if (is.matrix(data_input))
	{	data <- list(data_input)
	}
	else
	{	data <- data_input
	}
	if (is.null(rownames(data[[1]])))
	{	stop("rownames required")
	}
	d2 <- do.call(cbind,data)
	result <- .Call("ivc_multiple",
		d2,
		length(data),
		dim(data[[1]])[1],
		dim(data[[1]])[2],
		rownames(data[[1]]),
		treshold,
		treshold_multi
		)
	return(toBoolNet(names(result),result))
}

InferViaCorrelation_transitions <- function(data_input,treshold=0,alt=0)
#	data is list of state transitions
#	as created by timeseries_to_state_transitions
#	alt==0 -> modified correlation	(maucher et al. 2011 theroem 1)
#	alt==1 -> standard pearson correlation
{	if (is.matrix(data_input))
	{	data <- list(data_input)
	}
	else
	{	data <- data_input
	}
	if (is.null(rownames(data[[1]])))
	{       stop("rownames required")
	}
	st <- timeseries_to_state_transitions(data)

	d2 <- do.call(cbind,st)

	result <- .Call("ivc_ts",
			d2,
			length(st),	# count transition
			length(st[[1]]),	# size of transistion = 2*count genes
			rownames(data[[1]]),
			treshold,
			as.integer(alt)		# choose alternate correlation method
			)
	return(toBoolNet(rownames(data[[1]]),result))
}

InferViaCorrelation_transitions_higher_order <- function(data,treshold=0,order=1)
{	if (is.matrix(data_input))
	{	data <- list(data_input)
	}
	else
	{	data <- data_input
	}
	if (is.null(rownames(data[[1]])))
	{       stop("rownames required")
	}
	st <- timeseries_to_state_transitions(data)
	d2 <- do.call(cbind,st)
	result <- .Call("ivc_transitions_order",
			d2,
			length(st),     # count transition
			length(st[[1]]),        # size of transistion = 2*count genes
			rownames(data[[1]]),
			treshold,
			as.integer(order)
			)
	return(toBoolNet(rownames(data[[1]]),result))
}

ivc_last_corr <- function()
# debuggin returns all correlation coefficents of last InferViaCorrelation_transitions call
{	return(.Call("ivc_tr_last_corr"))
}

# debuggin purpose
ivc_last_ordered <- function()
{	return(.Call("ivc_tr_last_ordered"))
}

ivc_last_corr_ordered <- function()
{	return(.Call("ivc_tr_last_corr_ordered"))
}

toBoolNet <- function(genes,infered)
# converts result from InferViaCorrelation to BoolNet
# 	genes : list of all gene names
#	data : result from InferViaCorrelation
{	interactions <- lapply(genes, function(i)
		{	input <- unlist(lapply(infered[[i]], function(j)
				{       grep(j,genes)
				}))
#			func <- rep(0,2^length(input))    # TODO have and convert function
func <- c("?")
			expression <- paste0(
				c(
				"<f(",
				paste0(genes[as.numeric(input)],collapse=','),
				"){",
				paste0(func,collapse=''),
				"}>"),
				collapse='')
			list(input=input, func=func, expression=expression)
		})
	names(interactions) <- genes
	fixed <- rep(-1,length(genes))
	names(fixed) <- genes

	net <- list(genes=genes, interactions=interactions, fixed=fixed)
	class(net) <- "BooleanNetwork"

	return(net)
}

diffBN <- function(bn_1, bn_2)
# print differences; returns number of differences
{ 	if (class(bn_1)!="BooleanNetwork" || class(bn_2)!="BooleanNetwork")
	{	stop("class \"BooleanNetwork\" required")
	}
	diffs <- 0
	for ( i in setdiff(bn_1$genes,bn_2$genes))
	{	diffs <- diffs + 1
		cat(i,"missing in bn_2\n")
	}
	for (i in setdiff(bn_2$genes,bn_1$genes))
	{	diffs <- diffs + 1
		cat(i,"missing in bn_1\n")
	}
	for (i in intersect(bn_1$genes,bn_2$genes))
	{	d1_2 <- setdiff(bn_1$interactions[[i]]$input,bn_2$interactions[[i]]$input)
		d2_1 <- setdiff(bn_2$interactions[[i]]$input,bn_1$interactions[[i]]$input)

		if (length(d1_2) || length(d2_1))
		{	cat(i," : ")
			cat(paste(bn_1$interactions[[i]]$input)," <-> ")
			cat(paste(bn_2$interactions[[i]]$input),"\n")
			diffs <- diffs + length(d1_2) + length(d2_1)
		}
		# TODO compare func
	}

	return (diffs)
}

cmpBN <- function(bn_1, bn_2)
# compares two bn by dependencies
# positive is existing dependence; false is not existing dependence
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity
# returns named vector of true positives, false positives, true negatives, false negatives
# and all the other stuff
{ 	if (class(bn_1)!="BooleanNetwork" || class(bn_2)!="BooleanNetwork")
	{	stop("class \"BooleanNetwork\" required")
	}
	ret <- c(0,0,0,0)
	names(ret) <- c("TP","FP","TN","FN")
	P <- 0	# positive
	N <- 0	# negative

	all <- 1:length(bn_1$genes)
	for (i in intersect(bn_1$genes,bn_2$genes))
	{	ret["TP"] <- ret["TP"] + length(intersect(bn_1$interactions[[i]]$input,bn_2$interactions[[i]]$input))
		ret["FP"] <- ret["FP"] + length(setdiff(bn_2$interactions[[i]]$input,bn_1$interactions[[i]]$input))
		miss1 <- setdiff(all,bn_1$interactions[[i]]$input)
		miss2 <- setdiff(all,bn_2$interactions[[i]]$input)
		ret["TN"] <- ret["TN"] + length(intersect(miss1,miss2))
		ret["FN"] <- ret["FN"] + length(setdiff(miss2,miss1))
		P <- P + length(bn_1$interactions[[i]]$input)
		N <- N + length(bn_1$genes)-length(bn_1$interactions[[i]]$input)
	}

	ret2 <- rep(0,10)
	names(ret2) <- c("TPR","FPR","TNR","FNR","PPV","NPV","FDR","ACC","F1","MCC")
	ret2["TPR"] <- ret["TP"]/P
	ret2["FPR"] <- ret["FP"]/N
	ret2["TNR"] <- ret["TN"]/N
	ret2["FNR"] <- ret["FN"] / ( ret["TP"] + ret["FN"] )
	ret2["PPV"] <- ret["TP"] / ( ret["TP"] + ret["FP"] )
	ret2["NPV"] <- ret["TN"] / ( ret["TN"] + ret["FN"] )
	ret2["FDR"] <- ret["FP"] / ( ret["TP"] + ret["FP"] )
	ret2["ACC"] <- ( ret["TP"] + ret["TN"] ) / ( ret["TP"] + ret["FP"] + ret["FN"] + ret["TN"] )
	ret2["F1"]  <- 2 * ret["TP"] / ( 2 *  ret["TP"] + ret["FP"] + ret["FN"] )
	ret2["MCC"] <- ( ret["TP"] * ret["TN"] - ret["FP"] * ret["FN"] )/
		sqrt( ( ret["TP"] + ret["FP"] ) * ( ret["TP"] + ret["FN"] )
			* ( ret["TN"] + ret["FP"] ) * ( ret["TN"] + ret["FN"] ))

	return (c(ret,ret2))
}

timeseries_to_state_transitions <- function (timeseries)
# converts multiple timeseries to unique list of state-transistions
# state-transistion : first half is input; second half is output
{	c_genes=dim(timeseries[[1]])[1]
	ret <- list()
	li <- 1
	lapply(timeseries,function(i)
		{	samples=dim(i)[2]
			for (j in 1:(samples-1))
			{
				ret[[li]] <<- c(matrix(unlist(i),nrow=c_genes)[,j],
						matrix(unlist(i),nrow=c_genes)[,j+1])
				li <<- 1 + li
			}
		})
	return(unique(ret))
}

get_roc <- function(bn_1, bn_2)
# get receiver operating characteristics
# tpr	true positive rate
# fpr	false positive rate
{ 	if (class(bn_1)!="BooleanNetwork" || class(bn_2)!="BooleanNetwork")
	{	stop("class \"BooleanNetwork\" required")
	}
	ret <- c(0,0)
	names(ret) <- c("tpr","fpr")
	max <- 0
	for (i in intersect(bn_1$genes,bn_2$genes))
	{	ret["tpr"] <- ret["tpr"] + length(intersect(bn_1$interactions[[i]]$input,bn_2$interactions[[i]]$input))
		ret["fpr"] <- ret["fpr"] + length(setdiff(bn_2$interactions[[i]]$input,bn_1$interactions[[i]]$input))
		max <- max + c(
			length(bn_1$interactions[[i]]$input),
			length(bn_1$genes)-length(bn_1$interactions[[i]]$input)
			)
	}
	return (c(ret/max))
}

reorder_bn <- function(x)
# reorder function; assumes involved genes are already ordered
{	if (class(x)!="BooleanNetwork")
	{	stop("class \"BooleanNetwork\" required")
	}
	if (is.unsorted(x[[1]]))
	{	stop("base not ordered")
	}
	for (i in 1:length(x[[1]]))
	{	if (is.unsorted(x[[2]][[i]][[1]]))
		{	cat("sorting needed ",i,"\n")
			cat(x[[2]][[i]][[1]],"\n")
		}
	}

}

#******************** walsh - hadamard - transformation stuff ********************

is_power_of_2 <- function(x)
# true if x is power of two
{	return (1==sum(as.numeric(intToBits(x))))
}

log_2 <- function(x)
# returns dual logarithm for integers
{	if (x==1)
	{	return(1)
	}
	ret <- 1
	p <- 2
	while (x>p)
	{	ret <- ret + 1
		p <- p * 2
	}
	return (ret)
}

to_power_of_2 <- function(x)
# return 2^x
{	return(2^x)
}

# walsh indices w(0),w(1),w(2,1),w(2,2),w(3,1)..w(3,4),w(4,1)..w(4,8), ...
# walsh order    0    1     2       3      4  ..   7     8   ..  15

walsh_index_to_order <- function(a,b=NA)
# return order (0..n) from walsh indices
# input indices a,b or vector(2) containing a,b
{	if (is.na(b))
	{	b <- a[2]
		a <- a[1]
	}
	if (a==0) return(0)
	if (a==1) return(1)
	return( to_power_of_2(a-1) + b - 1 )
}

order_to_walsh_index <- function(o)
# walsh indices to order
{	if (o==0)
	{	return(0)
	}
	if (o==1)
	{	return(1)
	}
	a <- log_2(o+1)
	b <- o-to_power_of_2(a)/2+1
	return(c(a,b))
}

walsh_function <- function(order)
# creates walsh function, values 1,-1
{	if (order==0)
	{	return(1);
	}
	if (order==1)
	{	return(c(1,-1))
	}
	t <- walsh_function(as.integer(order/2))
	ret <- t
	if (order %% 2 == 1)
	{	t <- t * -1
	}
	return(c(ret,t))
}

walsh_hadamard_transformation <- function (x)
# vector of numeric of size 2^n
# returns same size vector of numeric
{	if (class(x)!="numeric")
	{	stop("walsh_hadamard_transformation : requires numeric argument")
	}
	if (!is_power_of_2(length(x)))
	{	stop("walsh_hadamard_transformation : size of argument must be power of two")
	}
	return (.Call("wht_r",
			x,
			length(x)
		))
}

diff_attractors_numeric <- function(a1,a2)
# compares AttractorInfo returns 0 if exact match
# any diff increases that value
# TODO ? some scoring values for various kinds of difference
{	if (class(a1)!="AttractorInfo" || class(a2)!="AttractorInfo")
	{	stop("parameters as class AttractorInfo required")
	}
	result <- 0
	result <- result + abs(length(a1$stateInfo$genes)-length(a2$stateInfo$genes))
	result <- result + abs(length(unique(a1$stateInfo$attractorAssignment))
				- length(unique(a2$stateInfo$attractorAssignment)))
	ca1 <- length(unique(a1$stateInfo$attractorAssignment))
	ca2 <- length(unique(a2$stateInfo$attractorAssignment))
	result <- result + abs(ca1 - ca2 )
	for (i in 1:ca1)
	{       found <- FALSE
		for (j in 1:ca2)
		{	s1 <- getAttractorSequence(a1,i)
			s2 <- getAttractorSequence(a2,j)
			if ((nrow(s1)==nrow(s2)) && (length( which(s1!=s2) )==0))
			{	found <- TRUE
			}
		}
		if (!found)
		{	result <- result + 1
		}
	}

	return(result)
}

diff_attractors <- function(a1,a2)
# compares AttractorInfo
# list differences
{	if (class(a1)!="AttractorInfo" || class(a2)!="AttractorInfo")
	{	stop("parameters as class AttractorInfo required")
	}
	# count involved genes
	cg1 <- length(a1$stateInfo$genes)
	cg2 <- length(a2$stateInfo$genes)
	if ( (cg1-cg2) != 0)
	{	cat("differing number of genes",diff_count_genes,"\n")
	}
	# count of found attractos
	ca1 <- length(unique(a1$stateInfo$attractorAssignment))
	ca2 <- length(unique(a2$stateInfo$attractorAssignment))
	if ( (ca1-ca2) != 0)
	{	cat("count attractors ",ca1," : ",ca2,"\n")
	}
	# cmp attractors
	for (i in 1:ca1)
	{	found <- FALSE
		s1 <- getAttractorSequence(a1,i)
		for (j in 1:ca2)
		{	s2 <- getAttractorSequence(a2,j)
			if ((nrow(s1)==nrow(s2)) && (length( which(s1!=s2))==0 ))
			{	found <- TRUE
			}
		}
		if (!found)
		{	cat("wrong attractor\n")
		}
	}
}	

# various helper stuff

create_all_possible_dependencies <- function(a)
# create all possible dependencies from vector of gene names
# in type used by reconstructNetwork
# depencies are given in named list of integers (indices into gene names)
{	ret <- lapply(a,function(i)
	{	match(c(a),a)
	})
	names(ret) <- a
	return (ret)
}

subtract_dependencies <- function(a,b)
# subtracts dependencies in b from a
{	if (!isTRUE(all.equal(names(a),names(b))))
	{	#stop("subtract_dependencies : colnames of a and b differing")
		stop("subtract_dependencies : colnames of a and b differing")
	}
	ret <- list()
	for (i in names(a))
	{	ret[[i]] <- setdiff(a[[i]],(b[[i]]))
	}
	return (ret)
}
