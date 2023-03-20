fx_loci_amplification_rate = function(ampseq_object, threshold = .65, chr_lengths = c(640851,
                                                                                       947102,
                                                                                       1067971,
                                                                                       1200490,
                                                                                       1343557,
                                                                                       1418242,
                                                                                       1445207,
                                                                                       1472805,
                                                                                       1541735,
                                                                                       1687656,
                                                                                       2038340,
                                                                                       2271494,
                                                                                       2925236,
                                                                                       3291936)){
  
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
  
  amplification_rate_per_locus = ggplot()+
    geom_segment(data = data.frame(x = 0, y = 1:14, xend =chr_lengths, yend = 1:14),
                 aes(x=x, y=y, xend=xend, yend=yend), alpha = .5)+
    geom_point(data = data.frame(Position = ampseq_object$markers$pos,
                                 Amplification_rate = 
                                   loci_performance$loci_ampl_rate,
                                 Chr = as.integer(gsub("(Pf3D7_)|(_v3)", "", ampseq_object$markers$chromosome))),
               aes(x = Position, y = Chr, color = Amplification_rate),
               pch = "|", size = 5)+
    theme_bw()+
    labs(title = "Amplification rate by locus",
         x = "Chromosome position",
         y = "Chromosome",
         color = NULL)+
    theme(legend.position = c(.9,.4))+
    scale_y_continuous(breaks = 1:14)+
    scale_color_continuous(type = "viridis")
  
  ampseq_loci_abd_table_discarded_loci = ampseq_loci_abd_table[, !(colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]])]
  ampseq_loci_abd_table = ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]]]
  
  markers = ampseq_object[["markers"]]
  
  discarded_markers = markers[loci_performance[["loci_ampl_rate"]] <= threshold ,]
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
                                           loci_performance = loci_performance_complete,
                                           markers = discarded_markers)
  ampseq_object[["plots"]][["all_loci_amplification_rate"]] = all_loci_performance_plot
  ampseq_object[["plots"]][["amplification_rate_per_locus"]] = amplification_rate_per_locus
  return(ampseq_object)
    
}
