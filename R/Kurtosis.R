# Purpose: Estimate kurtosis
# Updated: 180924

#' Estimate Kurtosis
#' 
#' Estimates the kurtosis of a distribution using a vector of observations. 
#' 
#' @param y Observations.
#' @param k0 Reference kurtosis. 
#' @param sig Significance level.
#' @param simple If TRUE, returns estimated size and SE only.
#' 
#' @return Matrix containing the estimated kurtosis, its standard error, the
#'   lower and upper confidence bounds, and a p value assessing the null
#'   hypothesis that the kurtosis is \code{k0}.
#' 
#' @importFrom stats pnorm qnorm
#' @export
#' 
#' @examples
#' \dontrun{
#' # Normal sample
#' y = rnorm(n=1e3);
#' # Test kurtosis = 3
#' Kurtosis(y=y,k0=3);
#' }

Kurtosis = function(y,k0=3,sig=0.05,simple=F){
  if(k0<=0){stop("Reference kurtosis must exceed zero.")};
  # Sample size
  n = length(y);
  # Mean
  m1 = mean(y);
  # Higher moments
  for(k in 2:8){
    cm = mean((y-m1)^k);
    assign(paste0("m",k),cm);
  }
  # Log Kurtosis
  lK = log(m4)-2*log(m2);
  # Asymptotic variance
  v = (m8)/(m4^2)-(4*m6)/(m2*m4)-(8*m3*m5)/(m4^2)+(4*m4)/(m2^2)+(16*m3^2)/(m2*m4)+(16*m2*m3^2)/(m4^2)-1;
  # Standard error
  SE = sqrt(v/n);
  # If simple, output
  if(simple){
    # Convert to original scale
    K = exp(lK);
    SE = K*SE;
    # Output
    Out = matrix(c(K,SE),nrow=1);
    colnames(Out) = c("Kurtosis","SE");
    rownames(Out) = c(1);
  } else {
  # Otherwise, calculate CI and p
    # Critical value
    t = qnorm(p=sig/2,lower.tail=F);
    # CI
    L = exp(lK-t*SE);
    U = exp(lK+t*SE);
    # P-value
    z = abs(lK-log(k0))/SE;
    p = 2*pnorm(q=z,lower.tail=F);
    # Convert to original scale
    K = exp(lK);
    SE = K*SE;
    # Format
    Out = matrix(c(K,SE,L,U,p),nrow=1);
    colnames(Out) = c("Kurtosis","SE","L","U","p");
    rownames(Out) = c(1);
  };
  return(Out);
}