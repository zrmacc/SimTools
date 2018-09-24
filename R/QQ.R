# Purpose: Function to generate QQ plots from summary counts
# Updated: 180921

#' Distribution of Negative Log P
#' 
#' Count the number of p values exceeding each break point
#' on the \eqn{-log_{10}} scale.
#' 
#' @param p Vector of p values.
#' @param U Upper limit of discretization grid.
#' @param B Bins.
#' 
#' @return A matrix containing the break points and the count of p-values
#'   exceeding each break point.
#' 
#' @export

allocateP = function(p,U=NULL,B=1000){
  # Observations
  n = length(p);
  if(is.null(U)){
    m = ceiling(log10(n));
  };
  # Log transform
  q = -log10(p);
  # Bin points
  b = seq(from=0,to=U,length.out=B+1);
  # Counts
  Counts = rep(0,B+1);
  Counts[1] = n;
  for(i in 2:(B+1)){
    keep = (q>b[i]);
    Counts[i] = sum(keep);
    q = q[keep];
  }
  # Output
  Out = cbind("Breaks"=b,"Counts"=Counts);
  rownames(Out) = seq(1:(B+1));
  return(Out);
}

#' Generate QQ Frame
#'
#' Generates points for the QQ plot. Assumes that \eqn{\beta_{0}=0}, and the
#' corresponding count is the total number of observations. Confidence intervals
#' are constructed using the beta distribution.
#' 
#' @param breaks Breakpoints.
#' @param counts Number of test statistics exceeding the breakpoint. 
#' 
#' @return A data.frame containing observed and expected uniform quantiles. 

qqFrame = function(breaks,counts){
  # Observations
  n = counts[1];
  # Survival function 
  S = counts/n;
  # Restrict to points with
  # 1. Positive survival
  # 2. Jump locations
  keep = (S>0)&(!duplicated(S));
  b = breaks[keep];
  S = S[keep];
  # Transform survival
  x = -log10(S);
  # Output
  df = data.frame("x"=x,"y"=b);
  return(df);
}

#' Construct QQ Plot
#' 
#' Constructs a uniform quantile-quantile plot using a sequence of breakpoints,
#' and the number of test statistics exceeding each breakpoint.
#' 
#' @param breaks Breakpoints.
#' @param counts Number of test statistics exceeding the breakpoint. 
#' @param sig Significance level for CIs.
#' @param color Point color.
#' 
#' @return A ggplot. 
#' 
#' @import ggplot2
#' @importFrom stats qbeta
#' @export

QQplot = function(breaks,counts,sig=0.05,color="gray"){
  # Construct plotting frame
  df = qqFrame(breaks=breaks,counts=counts);
  # Confidence bands
  n = counts[1];
  xk = (n+1)*10^(-breaks);
  # Lower limits
  Lk = -log10(qbeta(p=(1-sig/2),xk,n-xk+1));
  Uk = -log10(qbeta(p=(sig/2),xk,n-xk+1));
  # CI Frame
  cis = data.frame("x"=breaks,"L"=Lk,"U"=Uk);
  cis = cis[cis$x<=1.01*max(df$x),];
  
  # Initialize for check();
  x = y = L = U = NULL;
  # Plot
  q = ggplot() + theme_bw();
  q = q + geom_abline(intercept=0,slope=1,color="gray",linetype="dashed");
  q = q + geom_ribbon(data=cis,aes(x=x,ymin=L,ymax=U),alpha=0.2,fill=color);
  q = q + geom_point(data=df,aes(x=x,y=y),alpha=0.8,color=color);
  q = q + xlab("Expected Quantile") + ylab("Observed Quantile");
  return(q);
}