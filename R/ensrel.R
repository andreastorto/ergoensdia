#' Compute Ensemble Reliability Budget
#'
#' This function computes the ensemble reliability budget from observations
#' and observation-equivalent ensemble data.
#'
#' @param obs Vector of observed values (minimum 10 observations), of dimension n
#' @param ens Matrix containing the ensemble data at observation location, of dimension n x p (p being the ensemble size, minimum )
#' @param obs_err Value (constant) or vector of observation errors (either scalar or vector of dimension n)
#' @param plot Logical value for plotting the ensemble reliability budget as a barplot 
#' @param conf.l Confidence level to test the residual
#' @return A list with the ensemble reliability budget terms, including the p-value and the significance of the t-test for the null hypothesis (zero residual). A plot will be also produced (if plot == TRUE).
#' @export
#' @examples
#' ensrel(rnorm(100,mean=25,sd=1),array(rnorm(600,mean=25.5,sd=1.0),dim=c(100,6)),0.2,conf.l=0.99)

ensrel<-function(obs,ens,obs_err,plot=TRUE,conf.l=0.90){

        if ( ! is.vector(obs) && ! is.double(obs) ) stop("Provide obs as a vector of observations")


   	if ( ! is.matrix(ens) && ! is.array(ens) ) 
      	stop("Provide ens as a matrix or array (number of obs x ensemble size")

	ne<-length(obs_err)
   	if ( ne == 1 ) {
		   obs_err<-rep(obs_err,length(obs))
		   ne<-length(obs_err)
   	} 

	# Check the size of the problem

	no<-length(obs)
	nx<-dim(ens)[1]
	ns<-dim(ens)[2]

	if( no <  10 ) stop("Too few observations (less than 10)")
	if( nx != no || ne != no ) stop(paste("Mismatch in observation vector size",no,nx,ne))
	if( ns < 4   ) stop("Too few ensemble members (less than 4)")

	# Compute ensemble mean in observation space
	ensm<-rowMeans(ens)

	# Compute ensemble anomalies
	ensa<-array(NA,dim=dim(ens))
	for (j in 1:no) ensa[j,]<- ens[j,]-ensm[j]

	# Compute the budget terms
	n<-no
	n1<-n-1
	m<-ns
	m1<-m-1

	# Square departures
	depar2<-mean( (ensm-obs)^2 ) * (n/n1)

	# Square Bias
	bias2 <- ( mean( ensm-obs) )^2 * (n^2/(n*n1))

	# Ensemble variance
	ensvar <- sum( ensa^2 ) * (m+1)/(m*n*m1)

	# Observational error variance
	obserr2 <- mean(obs_err^2)

	# Squared residuals
	resid2 <- depar2-bias2-ensvar-obserr2

	# Test the significance of residuals
	ensas<-array(0,dim=length(ensm))
	for (j in 1:no) { ensas[j]<-sum(ensa[j,]^2) * (m+1)/(m*m1) }
	res2_l<-(ensm-obs)^2*(n/n1) - bias2 - ensas - obs_err^2
	ar<-t.test( res2_l, conf.level = conf.l)
	sign<-100*(1-ar$p.value)


	ensrel_out<-list( depar = depar2,
			  bias = bias2,
			  ensvar = ensvar,
			  obser = obserr2,
			  residual = resid2,
			  conf.int = ar$conf.int[1:2],
			  significance = sign )

	if ( plot ) {

		nams<-c("Square\nDeparture","Bias","Ensemble\nVariance",
			"Observation\nError","Residual")
		par(las=1,cex.main=1.2)
		vv<-c( depar2, bias2, ensvar, obserr2, resid2)
		aa<-range(vv[1:4],ar$conf.int)
		if(aa[1]>0) aa[1]<-0
		if(aa[2]<0) aa[2]<-0
		a<-barplot( vv, names.arg = nams, 
			ylim = aa,
			main="Ensemble reliability budget",
			col = c("lightblue", "mistyrose", "lightcyan",
                     "lavender", "cornsilk") )
		rr<-a[5]
		lines(c(rr,rr),ar$conf.int[1:2],lwd=2)
		for (i in 1:2) {
		  lines(c(rr-0.2,rr+0.2),c(ar$conf.int[i],ar$conf.int[i]),lwd=2)
		}

	}

	return (ensrel_out)
}
