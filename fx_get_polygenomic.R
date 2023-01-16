
fx_get_polygenomic = function(ampseq_object, strata, update_popsummary = T, na.rm = FALSE, filters = c('Buenaventura', 'Guapi', 'Quibdo')){
  
  library(Hmisc)
  
  ampseq_loci_abd_table = ampseq_object[["loci_abd_table"]]
  metadata = ampseq_object[["metadata"]]
  loci_performance = ampseq_object[["loci_performance"]]
  
  if(na.rm){
    ampseq_loci_abd_table = ampseq_loci_abd_table[!is.na(metadata[[strata]]) & !grepl('NA',metadata[[strata]]),]
    metadata = metadata[!is.na(metadata[[strata]]) & !grepl('NA',metadata[[strata]]),]
  }else{
    metadata[is.na(metadata[[strata]]) | !grepl('NA',metadata[[strata]]),][[strata]] = 'missing data'
  }
  
  if(!is.null(filters)){
    metadata = metadata[metadata[[strata]] %in% filters,]
    ampseq_loci_abd_table = ampseq_loci_abd_table[metadata[[strata]] %in% filters,]
  }
  
  polygenomic = NULL
  
  for(sample in rownames(ampseq_loci_abd_table)){
    polygenomic = rbind(polygenomic, data.frame(NPolyLoci = sum(grepl("_",ampseq_loci_abd_table[sample, ])),
               Polyloci = paste(names(ampseq_loci_abd_table[sample, ])[which(grepl("_",ampseq_loci_abd_table[sample, ]))], collapse = "/"),
               alleles_at_loci = paste(ampseq_loci_abd_table[sample, ][which(grepl("_",ampseq_loci_abd_table[sample, ]))], collapse = "/"),
               nalleles_per_loci = paste(
                 sapply(
                   sapply(
                     ampseq_loci_abd_table[sample, ][which(grepl("_", ampseq_loci_abd_table[sample, ]))], function(x){strsplit(x, "_")}),
                   function(x) length(x)),
                 collapse = "/"),
               coi = ifelse(sum(grepl("_",ampseq_loci_abd_table[sample, ])) == 0,
                            1,
                            max(unlist(sapply(
                              sapply(
                                ampseq_loci_abd_table[sample, ][which(grepl("_",ampseq_loci_abd_table[sample, ]))], function(x){strsplit(x, "_")}),
                              function(x) length(x)))))))
  }
  
  pop_summary = data.frame(
    pop = "Total",
    n = nrow(polygenomic),
    mean_coi = mean(polygenomic[["coi"]]),
    n_poly = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
    prop_poly = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                        n = nrow(polygenomic),
                        method = "exact")[1],
    prop_poly_lower = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                              n = nrow(polygenomic),
                              method = "exact")[2],
    prop_poly_upper = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                              n = nrow(polygenomic),
                              method = "exact")[3])
  
  
  for(pop in levels(as.factor(metadata[[strata]]))){
    pop_summary = rbind(pop_summary,
                        data.frame(
                          pop = pop,
                          n = nrow(polygenomic[metadata[[strata]] == pop,]),
                          mean_coi = mean(polygenomic[metadata[[strata]] == pop,][["coi"]]),
                          n_poly = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                          prop_poly = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                              n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                              method = "exact")[1],
                          prop_poly_lower = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                                    n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                    method = "exact")[2],
                          prop_poly_upper = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                                    n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                    method = "exact")[3]))
  }
  
  # Defining the proportion of polygenomic infections detected by loci --
  
  loci_performance[["prop_poly_detected"]] = apply(ampseq_loci_abd_table, 2, function(x) sum(grepl("_",x))/length(x))
  
  if(update_popsummary){
    ampseq_object[["metadata"]] = cbind(ampseq_object[["metadata"]], polygenomic)
    ampseq_object[["loci_performance"]] = loci_performance
    ampseq_object[["pop_summary"]] = if(is.null(ampseq_object[["pop_summary"]])){pop_summary}else{cbind(ampseq_object[["pop_summary"]], pop_summary)}
    return(ampseq_object)
  }else{
    return(pop_summary)
  }
  
}

