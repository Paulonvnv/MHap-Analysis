fx_sample_amplification_rate = function(ampseq_object, threshold = .8){
  
  metadata = ampseq_object[["metadata"]]
  ampseq_loci_abd_table = ampseq_object[["loci_abd_table"]]
  loci_performance = ampseq_object[["loci_performance"]]
  
  metadata[["sample_ampl_rate"]] = apply(ampseq_loci_abd_table, 1, function(x) 1 - sum(is.na(x))/length(x))
  
  
  all_samples_performance_plot = metadata %>% ggplot(aes(sample_ampl_rate))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = paste(sum(metadata[["sample_ampl_rate"]] > threshold), "Samples have an amplification rate >",threshold),
         x = "Amplification rate",
         y = "# of samples")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  ampseq_loci_abd_table_discarded_samples = ampseq_loci_abd_table[!(rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["samples"]]),]
  ampseq_loci_abd_table = ampseq_loci_abd_table[rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["samples"]],]
  
  metadata_complete = metadata
  metadata = metadata[metadata[["sample_ampl_rate"]] > threshold ,]
  
  
  loci_performance[["loci_ampl_rate2"]] = apply(ampseq_loci_abd_table, 2, function(x) 1- sum(is.na(x))/length(x))
  
  
  loci_performance_plot = loci_performance %>% ggplot(aes(loci_ampl_rate2))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = "Loci amplification rate",
         x = "Amplification rate",
         y = "# of Loci")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  
  ampseq_object[["loci_abd_table"]] = ampseq_loci_abd_table
  ampseq_object[["metadata"]] = metadata
  ampseq_object[["loci_performance"]] = loci_performance
  
  ampseq_object[["discarded_samples"]] = list(loci_abd_table = ampseq_loci_abd_table_discarded_samples,
                                           metadata = metadata_complete)
  
  ampseq_object[["plots"]][["loci_amplification_rate"]] = loci_performance_plot
  ampseq_object[["plots"]][["samples_amplification_rate"]] = all_samples_performance_plot
  
  return(ampseq_object)
  
  
}
