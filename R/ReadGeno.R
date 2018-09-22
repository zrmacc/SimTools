# Purpose: Function to import genotypes from plink 
# Updated: 180918

#' Draw Genotypes from Plink
#' 
#' Reads a random block of $m$ snps from a plink data set.
#' 
#' @param m Number of loci to import. 
#' @param stem Plink stem BED, BIM, FAM files.
#' @param IID Individuals to select. 
#' 
#' @importFrom snpStats read.plink
#' @export 

DrawGeno = function(m=1,stem,IID=NULL){
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
  # Snp information
  Snps = G$map[,c(1,2)];
  colnames(Snps) = c("Chr","Snp");
  rownames(Snps) = seq(1:nrow(Snps));
  # Genotypes
  H = matrix(as.numeric(G$genotypes),ncol=m);
  # Order
  if(!is.null(IID)){
    key = match(IID,G$fam$member);
  } else {
    key = seq(1:nrow(H));
  }
  H = H[key,,drop=F];
  # Convert to additive coding
  H = (3-H);
  H[H==3] = NA;
  # Label
  rownames(H) = G$fam$member[key];
  colnames(H) = Snps$Snp;
  # Output
  Out = list("Snps"=Snps,"Geno"=H);
}


#' Read Genotype Block from Plink
#' 
#' Reads a segment of genotypes from a plink file. 
#' 
#' @param L Beginning of segment.
#' @param U End of segment 
#' @param stem Plink stem BED, BIM, FAM files.
#' @param IID Individuals to select. 
#' 
#' @importFrom snpStats read.plink
#' @export 

GenoSeg = function(L=1,U=NULL,stem,IID=NULL){
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
  if(is.null(U)){U=M};
  if(U>M){stop("More loci requested than are available.")};
  
  # Snp block
  B = seq(from=L,to=U);
  b = length(B);
  # Import
  G = read.plink(bed=BED,bim=BIM,fam=FAM,select.snps=B);
  # Snp information
  Snps = G$map[,c(1,2)];
  colnames(Snps) = c("Chr","Snp");
  rownames(Snps) = seq(1:nrow(Snps));
  # Genotypes
  H = matrix(as.numeric(G$genotypes),ncol=b);
  # Order
  if(!is.null(IID)){
    key = match(IID,G$fam$member);
  } else {
    key = seq(1:nrow(H));
  }
  H = H[key,,drop=F];
  # Convert to additive coding
  H = (3-H);
  H[H==3] = NA;
  # Label
  rownames(H) = G$fam$member[key];
  colnames(H) = Snps$Snp;
  # Output
  Out = list("Snps"=Snps,"Geno"=H);
}