
fx_plot_network = function(relatedness_matrix,  
                        threshold,
                        metadata,
                        sample_id,
                        group_by,
                        levels,
                        colors){
  
  if(sum(is.na(metadata[[group_by]])) > 0){
    metadata[is.na(metadata[[group_by]]),][[group_by]] = 'missing data'
  }
  
  pairwise_relatednes_l = data.frame(Yi = rownames(relatedness_matrix), relatedness_matrix)
  
  pairwise_relatednes_l %<>% pivot_longer(cols = all_of(names(pairwise_relatednes_l)[-1]),
                                          names_to = "Yj",
                                          values_to = "r")
  
  pairwise_relatednes_l = pairwise_relatednes_l[!is.na(pairwise_relatednes_l$r),]
  
  pairwise_relatednes_l %<>% filter(Yi != Yj)
  
  variable = 'r'
  cols = c("Yi", "Yj")
  
  library(igraph)
  
  edges = NULL
  
  selected_edges <- pairwise_relatednes_l[pairwise_relatednes_l[[variable]] >= threshold,][,cols]
  for(edge in 1:nrow(selected_edges)){
    edges = c(edges,unlist(selected_edges[edge,]))
  }
  
  isolates = unique(unlist(pairwise_relatednes_l[,cols]))[!(unique(unlist(pairwise_relatednes_l[,cols])) %in% unique(edges))]
  
  network_object = graph(edges = edges, isolates = isolates, directed=F)
  
  node_colors = c(unique(edges), isolates)
  
  names(node_colors) = node_colors
  
  for(level in 1:length(levels)){
    node_colors[node_colors %in% metadata[metadata[[group_by]] == levels[level],][[sample_id]]] = colors[level]
  }
  
  plot_network = plot.igraph(network_object,
              vertex.color = node_colors,
              vertex.size=3,
              vertex.label.cex=0,
              vertex.label.dist=0.5, 
             # vertex.label.degree=-pi/2,
              vertex.label.color="black",
              vertex.label=NA,
              asp=0,
              xlim=c(-1,1),
              ylim=c(-1,1))
  
  return(list(network_object = network_object,
              plot_network = plot_network))
  
}