fx_loci_amplification_rate = function(ampseq_object, threshold = .65){
  
  ampseq_loci_abd_table = ampseq_object[["loci_abd_table"]]
  
  loci_performance = data.frame(loci = colnames(ampseq_loci_abd_table),
                                loci_ampl_rate = apply(ampseq_loci_abd_table, 2, function(x) 1 - sum(is.na(x))/length(x)))
  
  
  all_loci_performance_plot = loci_performance %>% ggplot(aes(loci_ampl_rate))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = paste(sum(loci_performance[["loci_ampl_rate"]] > threshold), "Loci have an amplification rate >",threshold),
         x = "Amplification rate",
         y = "# of Loci")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  ampseq_loci_abd_table_discarded_loci = ampseq_loci_abd_table[, !(colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]])]
  ampseq_loci_abd_table = ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]]]
  
  markers = ampseq_object[["markers"]]
  markers = markers[loci_performance[["loci_ampl_rate"]] > threshold ,]
  
  markers[["distance"]] = Inf
  
  for(chromosome in levels(as.factor(markers[["chromosome"]]))){
     for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
         markers[markers[["chromosome"]] == chromosome,][amplicon, "distance"] = markers[markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - markers[markers[["chromosome"]] == chromosome,][amplicon, "pos"]
       }
    }
  
  loci_performance_complete = loci_performance
  loci_performance = loci_performance[loci_performance[["loci_ampl_rate"]] > threshold ,]

  ampseq_object[["loci_abd_table"]] = ampseq_loci_abd_table
  ampseq_object[["markers"]] = markers
  ampseq_object[["loci_performance"]] = loci_performance
  ampseq_object[["discarded_loci"]] = list(loci_abd_table = ampseq_loci_abd_table_discarded_loci,
                                           loci_performance = loci_performance_complete)
  ampseq_object[["plots"]][["all_loci_amplification_rate"]] = all_loci_performance_plot
  
  return(ampseq_object)
    
}
