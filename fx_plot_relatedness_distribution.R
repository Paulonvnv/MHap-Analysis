
fx_plot_relatedness_distribution = function(relatedness_matrix = pairwise_relatedness,
                                            metadata = ampseq$metadata,
                                            Population = 'Population',
                                            fill_color = c("firebrick3", "firebrick1", "dodgerblue1", "dodgerblue3", "gold1", "gold3"),
                                            type_pop_comparison = 'within'){## c('within', 'between', 'both')
  
  pairwise_relatedness_l = data.frame(Yi = rownames(relatedness_matrix), relatedness_matrix)
  
  pairwise_relatedness_l %<>% pivot_longer(cols = all_of(names(pairwise_relatedness_l)[-1]),
                                          names_to = "Yj",
                                          values_to = "r")
  
  pairwise_relatedness_l = pairwise_relatedness_l[!is.na(pairwise_relatedness_l$r),]
  
  pairwise_relatedness_l %<>% filter(Yi != Yj)
  
  pairwise_relatedness_l$Pop_comparison = apply(pairwise_relatedness_l, 1, function(x){
    
    ifelse(metadata[metadata[['samples']] == x['Yi'],][[Population]] == 
       metadata[metadata[['samples']] == x['Yj'],][[Population]],
       metadata[metadata[['samples']] == x['Yi'],][[Population]],
       paste(sort(c(metadata[metadata[['samples']] == x['Yi'],][[Population]],
                   metadata[metadata[['samples']] == x['Yj'],][[Population]])), collapse = "-"))
  })
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    is.na(Pop_comparison) | grepl('NA', Pop_comparison) ~ "missing data",
    !(is.na(Pop_comparison) | grepl('NA', Pop_comparison)) ~ Pop_comparison
  ))
  
  pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within",
    Pop_comparison == 'missing data' ~ 'missing data'
  ))
  
  pairwise_relatedness_l %<>% filter(Pop_comparison != 'missing data', Type_of_comparison != 'missing data')
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                            apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')})))
  pairwise_relatedness_l$Type_of_comparison = factor(pairwise_relatedness_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  
  if(type_pop_comparison == 'within'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Within')%>%
      ggplot(aes(x = r, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$r), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~Pop_comparison,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
  }else if(type_pop_comparison == 'between'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Between')%>%
      ggplot(aes(x = r, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$r), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~Pop_comparison,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
  }else if(type_pop_comparison == 'both'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      ggplot(aes(x = r, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$r), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(Type_of_comparison~
                   Pop_comparison,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
  }

  
  pairwise_relatedness_distribution = list(relatedness = pairwise_relatedness_l,
                                           plot = plot_pairwise_relatedness_distribution)
  
  return(pairwise_relatedness_distribution)
  
}