#' Get a field from NetCDF file
#'
#' This function reads a netcdf variable
#' for quick use
#'
#' @param f Filename
#' @param v Variable name
#' @param start Optional array, specifying the "start"
#' @param count Optional array, specifying the "count"
#' @return The field
#' @export
#' @examples
#' temp<-get_fld("myfile.nc","temperature")
get_fld<-function(f,v,start=NULL,count=NULL){
    nc<-ncdf4::nc_open(f)
    if(is.null(start)||is.null(count)){
    var<-ncdf4::ncvar_get(nc,v)
    } else {
    var<-ncdf4::ncvar_get(nc,v,start=start,count=count)
    }
    ncdf4::nc_close(nc)
    return(var)
}

#' Put a field on an (existing) NetCDF file
#'
#' This function puts a netcdf variable
#' for quick use
#'
#' @param f Filename
#' @param v Variable name
#' @param var Field
#' @export
#' @examples
#' put_fld("myfile.nc","temperature",temp)
put_fld<-function(f,v,var){
    nc<-ncdf4::nc_open(f,write=TRUE)
    ncdf4::ncvar_put(nc,v,var)
    ncdf4::nc_close(nc)
}

#' Get an attribute from NetCDF file
#'
#' This function reads a netcdf attribute
#' for quick use
#'
#' @param f Filename
#' @param v Variable name
#' @param a Attribute name
#' @return The attribute
#' @export
#' @examples
#' temp_mv<-get_att("myfile.nc","temperature","missing_value")
get_att<-function(f,v,a){
    nc<-ncdf4::nc_open(f)
    att<-ncdf4::ncatt_get(nc,v,attname=a)
    ncdf4::nc_close(nc)
    return(att$value)
}

#' Compute the RMSE
#'
#' Compute the Root Mean Square Error
#'
#' @param x Verified field
#' @param y Verifying field
#' @return RMSE Value
#' @export
#' @examples
#' rmse(rnorm(100,mean=3,sd=2),rnorm(100,mean=3.5,sd=2.5))
rmse<-function(y,x){
 return( sqrt( mean ( (y-x)^2, na.rm=TRUE ) ) )
}

#' Decompose timeseries into trends, seasonal and residuals
#'
#' Decompose timeseries into interannual trends, annual and semi-annual signals and residuals
#'
#' @param data Timeseries to decompose
#' @param times Time, expressed in fractional years. If omitted, data are considered monthly
#' @return A list with the decomposed signal
#' @export
#' @examples
#' components( sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01 )
components<-function(data,times=NULL){
    nt<-length(data)
    if(times==NULL) {
	    mm<-(1:nt)/12
    } else { mm<-times }

    x1<-cos( mm*pi/0.5)
    x2<-sin( mm*pi/0.5)
    x3<-cos( mm*pi/0.25)
    x4<-sin( mm*pi/0.25)
    dd<-data.frame( y=data-mean(data), x0=mm, x1=x1, x2=x2, x3=x3, x4=x4 )
    lm.res<-lm( y ~ x0 + x1 + x2 + x3 + x4, data = dd )
    myt<-aov(lm.res)
    extrc<-anova(myt)["Sum Sq"]
    cdd<-sqrt( extrc/(nt-1) )
    cdd<-cdd[,1]
    zz <- atan2( lm.res$coeff[4], lm.res$coeff[3])
    zz2 <- atan2( lm.res$coeff[6], lm.res$coeff[5])
    pp<-zz ; pp2<-zz2
    zz <- abs( lm.res$coeff[4]/sin( zz ) )
    zz2<- abs( lm.res$coeff[6]/sin( zz2 ) )
    czz <- atan2( cdd[3], cdd[2])
    czz2 <- atan2( cdd[5], cdd[4])
    czz <- abs( cdd[3]/sin( czz ) )
    czz2<- abs( cdd[5]/sin( czz2 ) )
    fit<-lm.res$coeff[1]+lm.res$coeff[2]*mm+lm.res$coeff[3]*x1+
    lm.res$coeff[4]*x2+lm.res$coeff[5]*x3+lm.res$coeff[6]*x4
    fit1<-lm.res$coeff[1]+lm.res$coeff[2]*mm
    fit2<-lm.res$coeff[3]*x1+ lm.res$coeff[4]*x2+lm.res$coeff[5]*x3+lm.res$coeff[6]*x4
    return( list(a_a=zz,sa_a=zz2,trnd=lm.res$coeff[2]*12,a_p=pp,sa_p=pp2,
    fitted=fit,fitted_lt=fit1,fitted_seas=fit2,res=lm.res,acc=c(czz,czz2,cdd[1]*12)) )
}

#' Compute monthly climatology
#'
#' Compute monthly climatology from inter-annual monthly data.
#'
#' @param data Timeseries
#' @param mstop Abort if data do not cover years evenly (default to TRUE)
#' @return The 12-month climatology
#' @export
#' @examples
#' monthlyclim(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
monthlyclim <- function(data,mstop=TRUE){
        nr<-length(data)
        if( nr %% 12 != 0 && mstop ) stop("Monthlyclim works only with year-completed data")
        mc<-array(0,dim=c(12))
        for (i in 1:12){ mc[i]<-mean(data[seq(i,nr,by=12)],na.rm=TRUE) }
        return(mc)
}


#' Remove the mean
#'
#' Remove the mean from data
#'
#' @param data Timeseries
#' @return The anomaly timeseries
#' @export
#' @examples
#' demean(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
demean  <- function(data){
        return(data-mean(data))
}



#' Remove Trend
#'
#' Remove trend from data
#'
#' @param data Timeseries
#' @return The detrended timeseries
#' @export
#' @examples
#' rm_trend(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
rm_trend <- function(data){
        a<-components(data)
        newdata<-data-a$fitted_lt
        return(newdata)
}


#' Remove Seasonal signals
#'
#' Remove annual and semi-annual signals from data
#'
#' @param data Timeseries
#' @param useclim Use monthly climatology instead of fitting to estimate the seasonal signal (Default to FALSE)
#' @return The timeseries with seasonal signal removed
#' @export
#' @examples
#' rm_seas(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
rm_seas <- function(data,useclim=FALSE){
	if ( useclim ) {
		a<-monthlyclim(data)
        	a<-demean(a)
        	ny<-length(data)/12
        	newdata<-data-rep(a,ny)
        	return(newdata)
	} else {
        	a<-components(data)
        	newdata<-data-a$fitted_seas
        	return(newdata)
	}
}


#' Calculate periodogram
#'
#' Calculate (and plot) periodogram
#'
#' @param data Timeseries
#' @param plot Plot the frequency space data (default to TRUE)
#' @return The frequencies
#' @export
#' @examples
#' periodogram(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
periodogram <- function(data,plot=TRUE){
n<-length(data)
FF <- abs(fft(data)/sqrt(n))^2
nb<-(n/2)+1
nc<-(n/2)
P <- (4/n)*FF[1:nb] # Only need the first (n/2)+1 values of the FFT result.
f <- (0:nc)/n # this creates harmonic frequencies from 0 to .5
if(plot) plot(f*n, P, type="l",xlab="Months",log="y") #
return(list(x=f*n,y=P))
}

#' Smooth a timeseries
#'
#' Calculate a smoothed version of the timeseries using splines
#'
#' @param x time data
#' @param y timeseries to smooth
#' @param sp Spar parameter of splines (default to 0.2)
#' @return The smoothed timeseries
#' @export
#' @examples
#' smooth(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
smooth<-function(x,y,sp=0.2) {
x<-x[is.finite(y)]
y<-y[is.finite(y)]
if( (x[1]-x[2])> 0 ){
 s<-smooth.spline(x,y,spar=sp) ; return(predict(s,seq(x[1],x[length(x)],by=(x[1]-x[2])/3)))
} else {
 s<-smooth.spline(x,y,spar=sp) ; return(predict(s,seq(x[1],x[length(x)],by=(x[2]-x[1])/3)))
}
}

#' Lowpass or Highpass a timeseries
#'
#' Lowpass or Highpass a timeseries using a Butter filter
#'
#' @param data Timeseries
#' @param y timeseries to smooth
#' @param nm The cutoff frequency in units of time (i.e. use 12 with monthly data for 1-year frequency cutoff)
#' @return The filtered timeseries
#' @export
#' @examples
#' bf_lowpass (sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
#' bf_highpass(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
bf_lowpass<-function(data,nm){
  bf <- signal::butter(2, 1/nm, type="low")
  return(signal::filter(bf,data))
}
bf_highpass<-function(data,nm){
  bf <- signal::butter(2, 1/nm, type="high")
  return(signal::filter(bf,data))
}

#' Dates increment
#'
#' Increment dates by a certain number of days
#'
#' @param ymd Date in YYYYMMDD format (either string or numeric)
#' @param incr Increment in days
#' @return The new date
#' @export
#' @examples
#' dateincr(20120110,10)
dateincr<-function(ymd,incr=1){
   lint<-FALSE
   if(is.numeric(ymd)) {
        ymd<-sprintf("%d",ymd)
        lint<-TRUE
   }
   new<-format( as.Date(ymd,format= "%Y%m%d") + as.difftime(incr, unit="days") , "%Y%m%d" )
   if(lint) new<-as.integer(new)
   return(new)
}

#' Dates difference
#'
#' Calculate difference between two dates
#'
#' @param ymd1 Date in YYYYMMDD format (either string or numeric)
#' @param ymd2 Date in YYYYMMDD format (either string or numeric)
#' @return The days between ymd1 and ymd2
#' @export
#' @examples
#' datediff(20120120,20120110)
datediff<-function(ymd1,ymd2){
   if(is.numeric(ymd1)) {
        ymd1<-sprintf("%d",ymd1)
   }
   if(is.numeric(ymd2)) {
        ymd2<-sprintf("%d",ymd2)
   }
   new<-as.Date(ymd1,format= "%Y%m%d") - as.Date(ymd2,format= "%Y%m%d")
   return(as.integer(new))
}

#' Returns days of each month
#'
#' Returns days of each month, as a function of year
#'
#' @param year Optional (default non-leapfrog calendar)
#' @return The vector of days for the 12 months
#' @export
#' @examples
#' days_per_month(2005)
days_per_month<-function(year=NULL){
  mm<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  if(!is.null(year) && year%%4 == 0) mm[2]<-29
  return(mm)
}
