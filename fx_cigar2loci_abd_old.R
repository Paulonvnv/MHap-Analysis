
fx_cigar2loci_abd = function(cigar_object, min_abd = 10, min_ratio = .1, markers = NULL){
  
  cigar_table = cigar_object[["cigar_table"]]
  
  if(is.null(markers)){
    ampseq_loci_vector = unique(sapply(strsplit(colnames(cigar_table), ","), function(x) x[1]))
  }else{
    ampseq_loci_vector = markers[["amplicon"]]
  }

  ampseq_loci_abd_table = matrix(NA, nrow = nrow(cigar_table), ncol = length(ampseq_loci_vector), dimnames = list(rownames(cigar_table), ampseq_loci_vector))
  
  for(sample in rownames(ampseq_loci_abd_table)){
    for(locus in colnames(ampseq_loci_abd_table)){
      alleles = cigar_table[sample, grepl(paste("^",locus, sep = ""), colnames(cigar_table))]
      
      if(length(alleles) == 1){
        names(alleles) <- paste(locus, ",.", sep = "")
      }
      
      if(length(alleles[which(alleles > min_abd)]) == 1){
        ampseq_loci_abd_table[sample, locus] = paste(strsplit(names(alleles[which(alleles > min_abd)]), ",")[[1]][2], alleles[which(alleles > min_abd)], sep = ":")
      }else if(length(alleles[which(alleles > min_abd)]) > 1){
        alleles = sort(alleles, decreasing = T)
        ampseq_loci_abd_table[sample, locus] = gsub(paste(locus, ",", sep = ""), "", paste(paste(names(alleles[alleles/max(alleles) > min_ratio]), alleles[alleles/max(alleles) > min_ratio], sep = ":"), collapse = "_"))
        
      }
    }
    rm(list = c("sample", "locus", "alleles"))
  }
  
  rm(ampseq_loci_vector)
  
  # Removing controls
  
  controls_ampseq_loci_abd_table = ampseq_loci_abd_table[cigar_object[["metadata"]][cigar_object[["metadata"]][["typeofSamp"]] == "Controls",][["samples"]],]
  controls_metadata = cigar_object[["metadata"]][cigar_object[["metadata"]][["typeofSamp"]] == "Controls",]
  
  ampseq_loci_abd_table = ampseq_loci_abd_table[cigar_object[["metadata"]][cigar_object[["metadata"]][["typeofSamp"]] == "Samples",][["samples"]],]
  metadata = cigar_object[["metadata"]][cigar_object[["metadata"]][["typeofSamp"]] == "Samples",]
  
  return(list(loci_abd_table = ampseq_loci_abd_table,
              metadata = metadata,
              controls = list(loci_abd_table = controls_ampseq_loci_abd_table,
                              metadata = controls_metadata),
              markers = markers,
              loci_performance = NULL,
              pop_summary = NULL))
  

}