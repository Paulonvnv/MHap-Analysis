
fx_locus_diversity = function(locus, variance = TRUE){
  
  alleles = levels(as.factor(unlist(strsplit(gsub(":[0-9]+", "", locus), "_"))))
  
  n.all = length(alleles)
  
  n = sum(!is.na(unlist(strsplit(locus, "_"))))
  
  freq = sort(sapply(alleles, function(allele){
    sum(grepl(paste("(^|_)", allele, ":", sep = ""), locus))/n
  }), decreasing = T)
  
  sp2 <- sum(freq^2)
  H <- n * (1 - sp2)/(n - 1)
  na.e <- 1/(1-H)
  if (variance) {
    sp3 <- sum(freq^3)
    var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
    Nei.var.H <- (2*(n-1)/n^3)*((3-2*n)*(sp2^2)+2*(n-2)*sp3+sp2)
    return(c(n.all = n.all, na.e = na.e, Hexp = H, var.Hexp = var.H, Nei.var.Hexp = Nei.var.H))
  }
  else return(c(n.all = n.all, na.e = na.e, Hexp = H))
}
