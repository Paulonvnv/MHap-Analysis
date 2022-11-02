# Get allele frequencies--

fx_ampseq2loci = function(ampseq_object){

  ampseq_loci_abd_table = ampseq_object[["loci_abd_table"]]
  
  ampseq_allele_freq= NULL
  
  for(locus in colnames(ampseq_loci_abd_table)){
    alleles = levels(as.factor(unlist(strsplit(gsub(":[0-9]+", "", ampseq_loci_abd_table[,locus]), "_"))))
    ampseq_allele_freq[[locus]] = sort(sapply(alleles, function(allele){
      sum(grepl(paste("(^|_)", allele, ":", sep = ""), ampseq_loci_abd_table[,locus]))/sum(!is.na(unlist(strsplit(ampseq_loci_abd_table[,locus], "_"))))
    }), decreasing = T)
  }
  
  freq_table = matrix(NA, nrow = length(ampseq_allele_freq), ncol = max(sapply(ampseq_allele_freq, length)),
                      dimnames = list(names(ampseq_allele_freq),
                                      paste("Allele",seq(1,max(sapply(ampseq_allele_freq, length)),1), sep = "_")))
  
  for (locus in names(ampseq_allele_freq)){
    for(allele in 1:length(ampseq_allele_freq[[locus]])){
      freq_table[locus, allele] = ampseq_allele_freq[[locus]][allele]
    }
  }
  
  freq_table[is.na(freq_table)] = 0
  
  loci_table = gsub(":[0-9]+","" ,ampseq_loci_abd_table)
  
  for(locus in colnames(loci_table)){
    for (sample in rownames(loci_table)) {
      loci_table[sample,locus] = paste(which(names(ampseq_allele_freq[[locus]]) %in% strsplit(loci_table[sample,locus], "_")[[1]]), collapse = "_")
    }
  }
  
  loci_table[loci_table == ""] = NA
  
  return(list(loci_table = loci_table,
              freq_table = freq_table,
              ampseq_alleles_per_locus  = ampseq_allele_freq,
              markers = ampseq_object[["markers"]]))
  
}

