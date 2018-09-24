# Purpose: Estimate skewness
# Updated: 180924

#' Estimate Skewness
#' 
#' Estimates the skewness of a distribution using a vector of observations. 
#' 
#' @param y Observations.
#' @param s0 Reference skewness. 
#' @param sig Significance level.
#' @param simple If TRUE, returns estimated size and SE only.
#' 
#' @return Matrix containing the estimated skewness, its standard error, the
#'   lower and upper confidence bounds, and a p value assessing the null
#'   hypothesis that the skewness is \code{s0}.
#' 
#' @importFrom stats pnorm qnorm
#' @export
#' 
#' @examples
#' \dontrun{
#' # Exponential sample
#' y = rexp(n=1e3);
#' # Test skewness = 2
#' Skew(y=y,s0=2);
#' }

Skew = function(y,s0=0,sig=0.05,simple=F){
  # Sample size
  n = length(y);
  # Mean
  m1 = mean(y);
  # Higher moments
  for(k in 2:6){
    cm = mean((y-m1)^k);
    assign(paste0("m",k),cm);
  }
  # Skewness
  S = m3/(m2^(3/2));
  # Asymptotic variance
  v = m6/(m2^3)-(3*m3*m5)/(m2^4)+(9*m4*m3^2)/(m2^5)-(6*m4)/(m2^2)+(35*m3^2)/(4*m2^3)+9;
  # Standard error
  SE = sqrt(v/n);
  # If simple, output
  if(simple){
    # Format
    Out = matrix(c(S,SE),nrow=1);
    colnames(Out) = c("Skew","SE");
    rownames(Out) = c(1);
  } else {
  # Otherwise, calculate CI and p
    # Critical value
    t = qnorm(p=sig/2,lower.tail=F);
    # CI
    L = S-t*SE;
    U = S+t*SE;
    # P-value
    z = abs(S-s0)/SE;
    p = 2*pnorm(q=z,lower.tail=F);
    # Format
    Out = matrix(c(S,SE,L,U,p),nrow=1);
    colnames(Out) = c("Skew","SE","L","U","p");
    rownames(Out) = c(1);
  };
  # Output
  return(Out);
}