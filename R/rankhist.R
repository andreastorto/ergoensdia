#' Compute Rank Histogram
#'
#' This function computes the ensemble rank histogram budget from observations
#' and observation-equivalent ensemble data.
#'
#' @param obs Vector of observed values (minimum 10 observations), of dimension n
#' @param ens Matrix containing the ensemble data at observation location, of dimension n x p (p being the ensemble size, minimum )
#' @param plot Logical value for plotting the ensemble reliability budget as a barplot 
#' @param conf.l Confidence level to test the residual
#' @return A list with the ranks, including the RMSE of this histogram (w.r.t. to te perfect-ensemble value). A plot will be also produced (if plot == TRUE).
#' @export
#' @examples
#' rankhist(rnorm(100,mean=25,sd=1),array(rnorm(600,mean=25.5,sd=1.0),dim=c(100,6)))

rankhist<-function(obs,ens,plot=TRUE,conf.l=0.90){

        if ( ! is.vector(obs) && ! is.double(obs) ) stop("Provide obs as a vector of observations")

   	if ( ! is.matrix(ens) && ! is.array(ens) ) 
      	stop("Provide ens as a matrix or array (number of obs x ensemble size")

	# Check the size of the problem

	no<-length(obs)
	nx<-dim(ens)[1]
	ns<-dim(ens)[2]

	if( no <  10 ) stop("Too few observations (less than 10)")
	if( nx != no ) stop(paste("Mismatch in observation vector size",no,nx))
	if( ns < 4   ) stop("Too few ensemble members (less than 4)")

	# Compute the ranks and the perfect-ensemble value
	ranks <- apply(cbind(obs, ens), 1, rank, ties.method="random")[1, ]
	rankh<-array(0,dim=c(ns+1))
	for (i in 1:(ns+1)){ rankh[i]<-sum( ranks == i ) }
	print(rankh)
	perf<-no/(ns+1)

	# Compute the RMSE
	rmse<-sqrt ( mean ( (perf-ranks)^2 ) )

	rank_out<-list(   rankh = rankh,
			  rmse  = rmse,
			  nrnk  = ns+1 )

	if ( plot ) {

		par(las=1,cex.main=1.2)
		a<-barplot( rankh, names.arg = NULL,
			main="",border=NA,
			col = "white", xaxt="n",yaxt="n" )

		lines(c(-10,max(a)),c(perf,perf),lty=2,lwd=2)
		barplot( rankh, names.arg = sprintf("%02d",1:(ns+1)),
			main="Rank Histogram",
			col = "lightcyan", add = TRUE )
		legend("top",c(paste("RMSE =",round(rmse,1))),lty=0,bty="n")

	}

	return (rank_out)
}
