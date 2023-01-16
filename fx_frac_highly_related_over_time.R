
fx_plot_frac_highly_related_over_time = function(relatedness_matrix = pairwise_relatedness,
                                       metadata = ampseq$metadata,
                                       Population = c('Population', 'quarter_of_collection'),
                                       fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                       threshold = 0.99,
                                       type_pop_comparison = 'within'){## c('within', 'between', 'both')
  
  pairwise_relatedness_l = data.frame(Yi = rownames(relatedness_matrix), relatedness_matrix)
  
  pairwise_relatedness_l %<>% pivot_longer(cols = all_of(names(pairwise_relatedness_l)[-1]),
                                          names_to = "Yj",
                                          values_to = "r")
  
  pairwise_relatedness_l = pairwise_relatedness_l[!is.na(pairwise_relatedness_l$r),]
  
  pairwise_relatedness_l %<>% filter(Yi != Yj)
  
  # Giving wraning message, otherwise it is faster
  # pairwise_relatedness_l %<>% mutate(
  #   Population_Yi = metadata[grep(Yi, metadata[['samples']]),][[Population[1]]],
  #   Population_Yj = metadata[grep(Yj, metadata[['samples']]),][[Population[1]]],
  #   Date_Yi = metadata[grep(Yi, metadata[['samples']]),][[Population[2]]],
  #   Date_Yj = metadata[grep(Yj, metadata[['samples']]),][[Population[2]]])

  
  pairwise_relatedness_l[['Population_Yi']] = apply(pairwise_relatedness_l, 1, function(x){
    metadata[metadata[['samples']] == x['Yi'],][[Population[1]]]
  })
  
  pairwise_relatedness_l[['Population_Yj']] = apply(pairwise_relatedness_l, 1, function(x){
    metadata[metadata[['samples']] == x['Yj'],][[Population[1]]]
  })
  
  pairwise_relatedness_l[['Date_Yi']] = apply(pairwise_relatedness_l, 1, function(x){
    metadata[metadata[['samples']] == x['Yi'],][[Population[2]]]
  })
  
  pairwise_relatedness_l[['Date_Yj']] = apply(pairwise_relatedness_l, 1, function(x){
    metadata[metadata[['samples']] == x['Yj'],][[Population[2]]]
  })
  
  pairwise_relatedness_l %<>% filter(!is.na(Population_Yi),
                                     !is.na(Population_Yj),
                                     !is.na(Date_Yi),
                                     !grepl('NA',Date_Yi),
                                     !is.na(Date_Yj),
                                     !grepl('NA',Date_Yj))
  
  pairwise_relatedness_l %<>% mutate(Pop_Date_Yi = paste(Population_Yi, Date_Yi, sep = "_"),
                                    Pop_Date_Yj = paste(Population_Yj, Date_Yj, sep = "_"))
  
  
  plot_IBD_correlation_matrix = pairwise_relatedness_l %>%
    group_by(Pop_Date_Yi, Pop_Date_Yj) %>% 
    dplyr::summarise(freq = sum(r>=threshold),
                     n = n()) %>%
    group_by(Pop_Date_Yi, Pop_Date_Yj) %>% 
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
                           method = "exact")[3])%>%
    ggplot(aes(x = Pop_Date_Yi, y = Pop_Date_Yj, fill = prop)) + 
    geom_tile(alpha = .85)+
    scale_fill_gradient2(low = "white", high = "red", limit = c(0,.5), space = "Lab",
                         name=paste0('IBD >= ', threshold))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 12))
  

  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Population_Yi != Population_Yj ~ paste(Population_Yi, Population_Yj, sep = '_'),
    Population_Yi == Population_Yj ~ Population_Yi),
    Type_Pop_comparison = case_when(
    Population_Yi != Population_Yj ~ "Between",
    Population_Yi == Population_Yj ~ "Within"),
    Date_comparison = case_when(
      Date_Yi != Date_Yj ~ "Between",
      Date_Yi == Date_Yj ~ "Within")
    )
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(pairwise_relatedness_l$Pop_comparison)[!grepl("_", unique(pairwise_relatedness_l$Pop_comparison))]),
                                                            sort(unique(pairwise_relatedness_l$Pop_comparison)[grepl("_", unique(pairwise_relatedness_l$Pop_comparison))])))
  
  pairwise_relatedness_l$Type_Pop_comparison = factor(pairwise_relatedness_l$Type_Pop_comparison, levels =
                                                       c('Within', 'Between'))

  pairwise_relatedness_l %<>% filter(Date_comparison == 'Within')
  
  frac_highly_related = pairwise_relatedness_l %>%
    group_by(Pop_comparison, Date_Yi, Type_Pop_comparison) %>% 
    dplyr::summarise(freq = sum(r>=threshold),
                     n = n()) %>% group_by(Pop_comparison, Date_Yi, Type_Pop_comparison)%>%
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
                           method = "exact")[3])%>%
    mutate(lower = case_when(
      prop == 0 ~ 0,
      prop != 0 ~ lower),
      upper = case_when(
        prop == 0 ~ 0,
        prop != 0 ~ upper))
  
  if(type_pop_comparison == 'both'){
    
    plot_frac_highly_related = frac_highly_related %>%
      ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(Type_Pop_comparison~Pop_comparison)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'within'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Within')%>%
      ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~Pop_comparison)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Between')%>%
    ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~Pop_comparison)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }
  
  frac_highly_related = list(frac_highly_related = frac_highly_related,
                             plot_frac_highly_related = plot_frac_highly_related,
                             plot_IBD_correlation_matrix = plot_IBD_correlation_matrix)
  
  return(frac_highly_related)
  
}