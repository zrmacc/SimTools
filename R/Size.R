# Purpose: Estimate rejection probabilities
# Updated: 180918

#' Estimate Rejection Probability
#' 
#' Estimates the size of a hypothesis testing using a vector of p-values.
#' 
#' @param p Vector of p-values.
#' @param a Alpha level at which to assess size.
#' @param sig Significance level for CI.
#' 
#' @importFrom stats pnorm qnorm
#' @export

Size = function(p,a=0.05,sig=0.05){
  # Input check
  if(!is.vector(p)){stop("A numeric vector is expected for p.")};
  if((min(p)<0)||(max(p)>1)){stop("A vector of p-values is expected for p.")};
  # Missingness
  Miss = sum(is.na(p));
  if(Miss>0){stop("p should contain no missing data.")};
  
  # Obs
  n = length(p);
  # Estimate size
  m = mean(p<=a);
  # Standard error
  vM = (1/n)*max(m*(1-m),a*(1-a));
  SE = sqrt(vM);
  # Critical value
  t = qnorm(p=sig/2,lower.tail=F);
  # CI
  L = m-t*SE;
  U = m+t*SE;
  # P-value
  z = abs(m-a)/SE;
  p = 2*pnorm(q=z,lower.tail=F);
  # Output
  Out = matrix(c(m,SE,L,U,p),nrow=1);
  colnames(Out) = c("Point","SE","L","U","p");
  rownames(Out) = c(1);
  return(Out);
}