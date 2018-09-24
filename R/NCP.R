# Purpose: Estimate Non-Centrality Parameter
# Updated: 180918

#' Estimate Non-Centrality Parameter
#' 
#' Estimates the non-centrality parameter of the chi-square distribution using a
#' vector of p-values.
#' 
#' @param p Vector of p-values.
#' @param df Degrees of freedom.
#' @param sig Significance level for CI.
#' @param simple If TRUE, returns estimated size and SE only.
#' 
#' @return Matrix containing the estimated NCP, its standard error, the
#'   lower and upper confidence bounds, and a p-value assessing the null
#'   hypothesis the NCP is zero.
#' 
#' @importFrom stats pnorm qchisq qnorm var
#' @export 
#' 
#' @examples
#' \dontrun{
#' # Uniform p-values
#' p = runif(n=1e3);
#' # Non-centrality parameter
#' Ncp(p=p,df=1);
#' }

Ncp = function(p,df=1,sig=0.05,simple=F){
  # Input check
  if(!is.vector(p)){stop("A numeric vector is expected for p.")};
  if((min(p)<0)||(max(p)>1)){stop("A vector of p-values is expected for p.")};
  # Missingness
  Miss = sum(is.na(p));
  if(Miss>0){stop("p should contain no missing data.")};
  
  # Obs
  n = length(p);
  # Convert to chi-square statistics
  x = qchisq(p=p,df=df,lower.tail=F);
  # Estimate ncp
  d = mean(x)-df;
  # Standard error
  vD = (1/n)*var(x);
  SE = sqrt(vD);
  # If simple, output
  if(simple){
    # Format
    Out = matrix(c(d,SE),nrow=1);
    colnames(Out) = c("NCP","SE");
    rownames(Out) = c(1);
  } else {
  # Otherwise, calculate CI and p
    # Critical value
    t = qnorm(p=sig/2,lower.tail=F);
    # CI
    L = d-t*SE;
    U = d+t*SE;
    # P-value
    z = abs(d)/SE;
    p = 2*pnorm(q=z,lower.tail=F);
    # Format
    Out = matrix(c(d,SE,L,U,p),nrow=1);
    colnames(Out) = c("NCP","SE","L","U","p");
    rownames(Out) = c(1);
  };
  # Output
  return(Out); 
}