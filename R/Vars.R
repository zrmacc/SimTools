# Purpose: Compare the variances of two groups
# Updated: 180918

#' Compare Variances
#' 
#' Given two measurements on a single group, or measurements on two different
#' groups, calculates the difference and ratio of variances.
#' 
#' @param y1 Realizations at time $t_{1}$.
#' @param y0 Realizations at time $t_{0}$.
#' @param indep Are the measurements independent?
#' @param sig Significance level for CI.
#'   
#' @importFrom stats pnorm qnorm

comp.vars = function(y1,y0,indep=F,sig=0.05){
  # Output structure
  Out = array(0,dim=c(2,5));
  colnames(Out) = c("Point","SE","L","U","p");
  rownames(Out) = c(1:2);
  
  # Obs
  n1 = length(y1);
  n0 = length(y0);
  # First moments
  m11 = mean(y1);
  m01 = mean(y0);
  # Second moments
  m12 = var(y1);
  m02 = var(y0);
  # Fourth moment
  m14 = mean((y1-m11)^4);
  m04 = mean((y0-m01)^4);
  # Cross moments
  if(indep){
    s22 = m02*m12;
  } else {
    s22 = mean((y1-m11)^2*(y0-m01)^2);
  }
  
  # Critical value
  t = qnorm(p=sig/2,lower.tail=F);
  
  ## Difference
  D = m12-m02;
  # Standard error
  vD = (1/n0)*(m04-m02^2)+(1/n1)*(m14-m12^2)-(1/n0+1/n1)*(s22-m02*m12);
  SE = sqrt(vD);
  # CI
  L = D-t*SE;
  U = D+t*SE;
  # P-value
  z = abs(D)/SE;
  p = 2*pnorm(q=z,lower.tail=F);
  # Restults
  Out[1,] = c(D,SE,L,U,p);
  
  ## Ratio
  R = (m12/m02);
  # Standard error (log scale)
  vR = (1/n0)*m04/(m02^2)+(1/n1)*m14/(m12^2)-(1/n0+1/n1)*s22/(m02*m12);
  SE = sqrt(vR);
  # CI
  L = R*exp(-t*SE);
  U = R*exp(t*SE);
  # P-value
  z = abs(log(R))/SE;
  p = 2*pnorm(q=z,lower.tail=F);
  # Standard error (original scale)
  SE = R*SE;
  Out[2,] = c(R,SE,L,U,p);
  # Formatting
  Out = data.frame(Out);
  Out = cbind("Contrast"=c("Difference","Ratio"),Out);
  # Output
  return(Out);
}