log_scale_histogram = function(data, var, binwidth, group_by, levels, x_label, fill_color, y_breaks, ncol,
                               na.rm = FALSE,
                               filters = NULL){
  
  data = data[, c(var, group_by)]
  
  names(data) = c("variable", "groups")
  
  if(na.rm){
    data = data[!is.na(data[['groups']]) & !grepl('NA', data[['groups']]),]
  }else if(length(data[is.na(data[['groups']]) | grepl('NA', data[['groups']]),][['groups']])>0){
    data[is.na(data[['groups']]) | grepl('NA', data[['groups']]),][['groups']] = 'missing data'
  }
  
  if(!is.null(filters)){
    data = data[data[['groups']] %in% filters,]
  }
  
  bins = seq(min(data[["variable"]]), max(data[["variable"]]) + binwidth, binwidth)
  
  bins = data.frame(start = bins[-length(bins)], end = bins[-1])

  count_table = NULL
  
  for(group in levels(as.factor(data[["groups"]]))){
    temp_data = data[data[["groups"]]==group,]
    
    counts = sapply(1:nrow(temp_data), function(x) {
      which(bins[,1] <= temp_data[x,][["variable"]] & bins[,2] > temp_data[x,][["variable"]] )
    })
    
    for(bin in 1:nrow(bins)){
      count_table = rbind(count_table,
            data.frame(group = group,
                 bin = bins[bin,1],
                 Count = length(counts[counts == bin])))
    }
  }
  
  
  plot_histogram = count_table %>% ggplot(aes(x = bin, y = log10(Count + 1), fill = factor(group, levels = levels)))+
    geom_col(alpha = .7)+
    facet_wrap(.~factor(group,
                        levels = levels), ncol = ncol)+
    scale_y_continuous(breaks = log10(y_breaks + 1), labels = y_breaks)+
    labs(y = "Count",
         x = x_label)+
    scale_fill_manual(values = fill_color)+
    theme_bw()+
    theme(legend.position = "none")
  
  return(plot_histogram)
  
}
