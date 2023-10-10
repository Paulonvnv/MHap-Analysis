
fx_plot_frac_highly_related = function(relatedness_matrix = pairwise_relatedness,
                                      metadata = ampseq$metadata,
                                      Population = 'Population',
                                      fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                      threshold = 0.99,
                                      type_pop_comparison = 'between',
                                      pop_levels = NULL){
  
  pairwise_relatedness_l = data.frame(Yi = rownames(relatedness_matrix), relatedness_matrix)
  
  pairwise_relatedness_l %<>% pivot_longer(cols = all_of(names(pairwise_relatedness_l)[-1]),
                                          names_to = "Yj",
                                          values_to = "r")
  
  pairwise_relatedness_l = pairwise_relatedness_l[!is.na(pairwise_relatedness_l$r),]
  
  pairwise_relatedness_l %<>% filter(Yi != Yj)
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('samples', Population)], by.x = 'Yi', by.y = 'samples', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-4], 'Yi_Population')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('samples', Population)], by.x = 'Yj', by.y = 'samples', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Yj_Population')
  
  pairwise_relatedness_l %<>% filter(!is.na(Yi_Population), !is.na(Yj_Population))
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Yi_Population == Yj_Population ~ Yi_Population,
    Yi_Population != Yj_Population ~ max(Yi_Population, Yj_Population)
  ))
  
  pairwise_relatedness_l$Pop_comparison = apply(pairwise_relatedness_l, 1, function(x){
    
    ifelse( x['Yi_Population'] == x['Yj_Population'],
            x['Yi_Population'],
            paste(sort(c(x['Yi_Population'], x['Yj_Population'])), collapse = "-"))
  })
  
  # if(sum(is.na(pairwise_relatedness_l$Pop_comparison)) + sum(grepl('NA', pairwise_relatedness_l$Pop_comparison)) >= 1){
  #   pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
  #     is.na(Pop_comparison) | grepl('NA', Pop_comparison) ~ "missing data",
  #     !(is.na(Pop_comparison) | grepl('NA', Pop_comparison)) ~ Pop_comparison
  #   ))
  #   
  #   pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
  #     grepl("-",Pop_comparison) ~ "Between",
  #     !grepl("-",Pop_comparison) ~ "Within",
  #     Pop_comparison == 'missing data' ~ 'missing data'
  #   ))
  #   
  #   pairwise_relatedness_l %<>% filter(Pop_comparison != 'missing data', Type_of_comparison != 'missing data')
  # }else{
  #   
  pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within"
  ))
  
  # }
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                            apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')})))
  pairwise_relatedness_l$Type_of_comparison = factor(pairwise_relatedness_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  highly_related_table = pairwise_relatedness_l %>%
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
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
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
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
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
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
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