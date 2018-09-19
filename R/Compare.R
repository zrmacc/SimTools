# Purpose: Main two-group comparison function
# Updated: 180918

#' Compare Cumulants
#' 
#' Given two measurements on a single group, or measurements on two different 
#' groups, compares the group means or variances.
#' 
#' @param y1 Measurements at time $t_{1}$.
#' @param y0 Measurements at time $t_{0}$.
#' @param compare Either means or vars.
#' @param indep Are the measurements independent?
#' @param sig Significance level for CI.
#' 
#' @export

compCumulants = function(y1,y0,compare="means",indep=F,sig=0.05){
  # Input check
  if(!is.vector(y1)){stop("A numeric vector is expected for y1.")};
  if(!is.vector(y0)){stop("A numeric vector is expected for y0.")};
  # Missingness
  Miss = sum(is.na(y1))+sum(is.na(y0));
  if(Miss>0){stop("y1 and y0 should contain no missing data.")};
  # Comparison
  Choices = c("means","vars");
  if(!(compare%in%Choices)){stop("Select comparison from among: means, vars, efficiency.")};
  # Dependence
  if((!indep)&&(length(y1)!=length(y0))){stop("If comparing two measurements on the same group, y1 and y0
                                              are expected to have the same length.")};
  # Means
  if(compare=="means"){
    Out = comp.means(y1=y1,y0=y0,indep=indep,sig=sig);
  } else if(compare=="vars"){
    Out = comp.vars(y1=y1,y0=y0,indep=indep,sig=sig);
  }
  return(Out);
}

########################
# Efficiency
########################

#' Compare Efficiency
#' 
#' Given p-values from two tests of the same hypothesis on the same data,
#' or on two different data sets, compares the efficiency using the 
#' difference and ratio of non-centrality parameters. 
#' 
#' @param p1 P-values from test 1.
#' @param p0 P-values from test 0.
#' @param indep Are the measurements independent?
#' @param df Degress of freedom, if comparing the efficiency of hypothesis
#'   tests.
#' @param sig Significance level for CI.
#' 
#' @importFrom stats qchisq
#' @export

compEfficiency = function(p1,p0,df=1,indep=F,sig=0.05){
  # Input check
  if(!is.vector(p1)){stop("A numeric vector is expected for p1.")};
  if(!is.vector(p0)){stop("A numeric vector is expected for p0.")};
  # Missingness
  Miss = sum(is.na(p1))+sum(is.na(p0));
  if(Miss>0){stop("p1 and p0 should contain no missing data.")};
  # Validity
  if((min(p1)<0)||(min(p0)<0)||(max(p1)>1)||(max(p0)>1)){
    stop("Vectors of p-values are expected for p1 and p0.")};
  # Dependence
  if((!indep)&&(length(p1)!=length(p0))){
    stop("If comparing two measurements on the same group, p1 and p0 are expected to have the same length.")};
  # Convert to chi-square statistics
  y1 = qchisq(p=p1,df=df,lower.tail=F)-df;
  y0 = qchisq(p=p0,df=df,lower.tail=F)-df;
  # Compare non-centrality parameters
  Out = comp.means(y1=y1,y0=y0,indep=indep,sig=sig);
  return(Out);
}
