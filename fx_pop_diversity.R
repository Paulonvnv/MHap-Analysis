fx_pop_diversity = function(ampseq_object, strata){
  
  source("fx_loci_diversity.R")
  
  ampseq_loci_abd_table = ampseq_object[["loci_abd_table"]]
  metadata = ampseq_object[["metadata"]]
  
  mlg = as.factor(apply(gsub(":[0-9]+", "",ampseq_loci_abd_table[metadata[["coi"]] == 1, ]), 1, function(sample) paste(sample, collapse = "_")))
  
  n_mlg = nlevels(mlg) # Richness: number of different variants (species, genus, families, ASVs)
  p = summary(mlg)/length(mlg)# Frequency of each variant
  sp2 = sum(p^2)
  D = length(mlg) * (1 - sp2)/(length(mlg) - 1) # Simpson: probability of choosing two different variants
  S.e = 1/(1-D) # effective richness: The number of variants that explain the diversity of Simpson
  H = -sum(p*log(p)) # Shannon
  E = H/log(n_mlg) # Evenness: from 0 to 1
  
  loci_diversity = fx_loci_diversity(ampseq_loci_abd_table, variance = F)
  
  pop_diversity = data.frame(pop = "Total",
                             Richness = n_mlg,
                             Effective.richness = S.e,
                             Simpson = D,
                             Shannon = H,
                             Evenness = E,
                             n.all = mean(loci_diversity$n.all),
                             na.e = mean(loci_diversity$na.e),
                             Hexp = mean(loci_diversity$Hexp)
                             )
  
  for(pop in levels(as.factor(metadata[[strata]]))){
    
    if(is.null(nrow(ampseq_loci_abd_table[metadata[[strata]] == pop & metadata[["coi"]] == 1, ]))){
      
      pop_diversity = rbind(pop_diversity,data.frame(pop = pop,
                                                     Richness = 1,
                                                     Effective.richness = 1,
                                                     Simpson = 0,
                                                     Shannon = -sum(1*log(1)),
                                                     Evenness = -sum(1*log(1))/log(1),
                                                     n.all = 1,
                                                     na.e = 1,
                                                     Hexp = 0))
    }else{
        mlg = as.factor(apply(gsub(":[0-9]+", "",ampseq_loci_abd_table[metadata[[strata]] == pop & metadata[["coi"]] == 1, ]), 1, function(sample) paste(sample, collapse = "_")))
        n_mlg = nlevels(mlg) # Richness: number of different variants (species, genus, families, ASVs)
        p = summary(mlg)/length(mlg)# Frequency of each variant
        sp2 = sum(p^2)
        D = length(mlg) * (1 - sp2)/(length(mlg) - 1) # Simpson: probability of choosing two different variants
        S.e = 1/(1-D) # effective richness: The number of variants that explain the diversity of Simpson
        H = -sum(p*log(p)) # Shannon
        E = H/log(n_mlg) # Evenness: from 0 to 1
        
        loci_diversity = fx_loci_diversity(ampseq_loci_abd_table[metadata[[strata]] == pop, ], variance = F)
        
        pop_diversity = rbind(pop_diversity,data.frame(pop = pop,
                                                       Richness = n_mlg,
                                                       Effective.richness = S.e,
                                                       Simpson = D,
                                                       Shannon = H,
                                                       Evenness = E,
                                                       n.all = mean(loci_diversity$n.all),
                                                       na.e = mean(loci_diversity$na.e),
                                                       Hexp = mean(loci_diversity$Hexp)))
    }
    
  }
  

  return(pop_diversity)
  
}
