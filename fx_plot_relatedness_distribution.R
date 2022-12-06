
fx_plot_relatedness_distribution = function(relatedness_matrix = pairwise_relatedness,
                                            metadata = ampseq$metadata,
                                            Population = 'Population',
                                            fill_color = c("firebrick3", "firebrick1", "dodgerblue1", "dodgerblue3", "gold1", "gold3")){
  
  pairwise_relatednes_l = data.frame(Yi = rownames(relatedness_matrix), relatedness_matrix)
  
  pairwise_relatednes_l %<>% pivot_longer(cols = all_of(names(pairwise_relatednes_l)[-1]),
                                          names_to = "Yj",
                                          values_to = "r")
  
  pairwise_relatednes_l = pairwise_relatednes_l[!is.na(pairwise_relatednes_l$r),]
  
  pairwise_relatednes_l %<>% filter(Yi != Yj)
  
  pairwise_relatednes_l$Pop_comparison = apply(pairwise_relatednes_l, 1, function(x){
    
    ifelse(metadata[metadata[['samples']] == x['Yi'],][[Population]] == 
       metadata[metadata[['samples']] == x['Yj'],][[Population]],
       metadata[metadata[['samples']] == x['Yi'],][[Population]],
       paste(sort(c(metadata[metadata[['samples']] == x['Yi'],][[Population]],
                   metadata[metadata[['samples']] == x['Yj'],][[Population]])), collapse = "-"))
  })
  
  
  pairwise_relatednes_l %<>% mutate(Type_of_comparisson = case_when(
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within"
  ))
  
  pairwise_relatednes_l$Pop_comparisson = factor(pairwise_relatednes_l$Pop_comparisson,
                                                 levels = c(sort(unique(metadata[[Population]])),
                                                            apply(combn(unique(metadata[[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')})))
  pairwise_relatednes_l$Type_of_comparisson = factor(pairwise_relatednes_l$Type_of_comparisson, levels =
                                                       c('Within', 'Between'))
  
  plot_pairwise_relatedness_distribution = pairwise_relatednes_l %>%
    ggplot(aes(x = r, fill = Pop_comparison)) +
    geom_histogram(position = "stack", alpha = .7)+
    geom_vline(xintercept = mean(pairwise_relatednes_l$r), linetype = 2)+
    scale_fill_manual(values = fill_color)+
    theme_bw()+
    facet_wrap(Type_of_comparisson~
                 Pop_comparisson,
               scales = "free_y")+
    labs(y = "Count",
         x = "Relatedness",
         fill = "Populations")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          legend.position = "none")
  
  return(plot_pairwise_relatedness_distribution)
  
}