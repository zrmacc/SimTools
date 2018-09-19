# Purpose: Function to import genotypes from plink 
# Updated: 180918

#' Read Genotypes from Plink
#' 
#' Reads a random block of $m$ snps from a plink data set.
#' 
#' @param m Number of loci to import. 
#' @param stem Plink stem BED, BIM, FAM files.
#' 
#' @importFrom snpStats read.plink
#' @export 

ReadGeno = function(m=1,stem){
  # Plink files
  BED = paste0(stem,".bed");
  BIM = paste0(stem,".bim");
  FAM = paste0(stem,".fam");
  # Existence
  if(!file.exists(BIM)){stop("BIM file DNE.")};
  # Available loci
  cmd = paste0("wc -l ",BIM);
  M = trimws(system(command=cmd,intern=T));
  # Parse
  M = as.integer(gsub(pattern="^([0-9]+)\\ .*",replacement="\\1",x=M));
  if(m>M){stop("More loci requested than are available.")};
  
  # Draw starting location
  Draw = sample(x=M-m,size=1);
  # Snp block
  B = seq(from=Draw,to=(Draw+m-1));

  # Import
  G = read.plink(bed=BED,bim=BIM,fam=FAM,select.snps=B);
  G = matrix(as.numeric(G$genotypes),ncol=m);
  # Convert to additive coding
  G = (3-G);
  G[G==3] = NA;
  return(G);
}