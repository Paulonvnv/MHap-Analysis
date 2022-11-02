fx_loci_diversity = function(loci_abd_table = NULL, variance = TRUE){
  
  source("fx_locus_diversity.R")
  
  loci_diversity = data.frame(t(sapply(colnames(loci_abd_table), function(locus) fx_locus_diversity(loci_abd_table[,locus], variance = variance))))

  return(loci_diversity)


}
