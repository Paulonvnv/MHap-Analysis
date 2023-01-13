
fx_plot_frac_highly_related = function(relatedness_matrix = pairwise_relatedness,
                                      metadata = ampseq$metadata,
                                      Population = 'Population',
                                      fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                      threshold = 0.99,
                                      type_pop_comparison = 'between'){
  
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
  
  
  pairwise_relatednes_l %<>% mutate(Type_of_comparison = case_when(
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within"
  ))
  
  pairwise_relatednes_l$Pop_comparison = factor(pairwise_relatednes_l$Pop_comparison,
                                                 levels = c(unique(metadata[[Population]]),
                                                            sort(apply(combn(unique(metadata[[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')}))))
  pairwise_relatednes_l$Type_of_comparison = factor(pairwise_relatednes_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  highly_related_table = pairwise_relatednes_l %>%
    group_by(Type_of_comparison, Pop_comparison) %>% 
    dplyr::summarise(freq = sum(r>=threshold),
                     n = n()) %>% group_by(Type_of_comparison, Pop_comparison)%>%
    mutate(prop = binconf(freq,
                          n,
                          alpha = 0.05,
                          method = "exact")[1],
           lower = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[2],
           upper = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[3]) %>%
    mutate(lower = case_when(
      prop == 0 ~ 0,
      prop != 0 ~ lower),
      upper = case_when(
        prop == 0 ~ 0,
        prop != 0 ~ upper))
  
  if(type_pop_comparison == 'within'){
    
    plot_frac_highly_related = highly_related_table %>%
      filter(Type_of_comparison == 'Within') %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = highly_related_table %>%
      filter(Type_of_comparison == 'Between') %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'both'){
    
    plot_frac_highly_related = highly_related_table %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }
  
  
  
  frac_highly_related = list(highly_related_table = highly_related_table,
                             plot = plot_frac_highly_related)
  
  return(frac_highly_related)
  
}