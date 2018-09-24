# Purpose: Characterize the cumulants of a sample
# Updated: 180924

#' Characterize Cumulants
#' 
#' Estimates the cumulants of a distribution using a vector of observations. 
#' 
#' @param y Observations.
#' @param sig Significance level.
#' @param report Display result? 
#' 
#' @return Data.frame containing the estimated mean, variance, skewness, and kurtosis.
#' 
#' @importFrom stats qnorm
#' @export
#' 
#' @examples
#' \dontrun{
#' # Normal sample
#' y = rnorm(n=1e3);
#' # Characterize
#' K = Characterize(y=y);
#' }

Characterize = function(y,sig=0.05,report=T){
  # Obs
  n = length(y);
  # Mean
  m1 = mean(y);
  # Variance
  m2 = var(y);
  # Fourth moment
  m4 = mean((y-m1)^4);
  # Standard errors
  v1 = m2;
  SE1 = sqrt(v1/n);
  v2 = m4/(m2^2)-1;
  SE2 = sqrt(v2/n);
  # Critical value
  t = qnorm(p=sig/2,lower.tail=F);
  # CI for mean
  L1 = m1-t*SE1;
  U1 = m1+t*SE1;
  # CI for var
  L2 = m2*exp(-t*SE2);
  U2 = m2*exp(+t*SE2);
  # Convert to original scale
  SE2 = m2*SE2;
  # Mean
  M = matrix(c(m1,SE1,L1,U1),nrow=1);
  # Variance
  V = matrix(c(m2,SE2,L2,U2),nrow=1);
  # Skewness
  S = Skew(y=y,sig=sig)[,1:4];
  # Kurtosis
  K = Kurtosis(y=y,sig=sig)[,1:4];
  # Format
  Out = data.frame(rbind(M,V,S,K));
  colnames(Out)[1] = "Point";
  rownames(Out) = seq(1:4);
  Out$Cumulant = c("Mean","Var","Skew","Kurt");
  Out = Out[,c(5,1:4)];
  # Report
  if(report){
    Disp = Out;
    aux = function(v){
      if(is.numeric(v)){return(signif(v,digits=3))}
      else{return(v)};
    };
    Disp[] = lapply(Disp,aux);
    print(Disp);
  }
  # Output
  return(Out);
}
