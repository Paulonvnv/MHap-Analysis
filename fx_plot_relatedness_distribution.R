
fx_plot_relatedness_distribution = function(relatedness_matrix = pairwise_relatedness,
                                            metadata = ampseq$metadata,
                                            Population = 'Population',
                                            fill_color = c("firebrick3", "firebrick1", "dodgerblue1", "dodgerblue3", "gold1", "gold3"),
                                            type_pop_comparison = 'within',
                                            ncol = 4,
                                            pop_levels = NULL
                                            ){## c('within', 'between', 'both')
  
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
  
  
  
  if(type_pop_comparison == 'within'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Within')%>%
      ggplot(aes(x = r, fill = factor(Pop_comparison, levels = pop_levels))) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$r), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol,
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
      ggplot(aes(x = r, fill = factor(Pop_comparison, levels = pop_levels))) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$r), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol,
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
