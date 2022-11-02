fx_pairwise_relatedness = function(loci_object){
  
  loci_table = loci_object[["loci_table"]]
  freq_table = loci_object[["freq_table"]]
  markers = loci_object[["markers"]]
  
  pairs = rbind(cbind(rownames(loci_table), rownames(loci_table)),
                t(combn(rownames(loci_table), 2)))

  relatedness_matrix = matrix(NA,
                              nrow = nrow(loci_table),
                              ncol = nrow(loci_table),
                              dimnames = list(rownames(loci_table),
                                              rownames(loci_table)))
  
  fx_get_relatedness = function(Yi, Yj, freq_table, markers){
    
    Ys = cbind(Yi, Yj)
    
    rownames(Ys) = colnames(loci_table)
    
    Ys =Ys[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    fs = freq_table[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers = markers[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers[["distance"]] = Inf
    
    for(chromosome in levels(as.factor(no_na_markers[["chromosome"]]))){
      for(amplicon in 1:(nrow(no_na_markers[no_na_markers[["chromosome"]] == chromosome,])-1)){
        no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "distance"] = no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "pos"]
      }
    }
    
    ds = no_na_markers$distance
    
    estimate = paneljudge::estimate_r_and_k(fs = fs,
                                            ds = ds,
                                            Ys = Ys,
                                            warn_fs = F)
    
    return(estimate["rhat"])
  }
  
  
  for(pair in 1:nrow(pairs)){
    
    if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
      Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
      Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
      
      estimate = NULL
      for(haplotype_i in 1:ncol(Yi)){
        for(haplotype_j in 1:ncol(Yj)){
          estimate = c(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj[,haplotype_j], freq_table, markers))
        }
      }
      
      relatedness_matrix[pairs[pair, 2],pairs[pair, 1]] = max(estimate)
      
    }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
      Yi = as.numeric(loci_table[pairs[pair, 1],])
      Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
      
      estimate = NULL
      for(haplotype_j in 1:ncol(Yj)){
        estimate = c(estimate, fx_get_relatedness(Yi, Yj[,haplotype_j], freq_table, markers))
      }
      
      relatedness_matrix[pairs[pair, 2],pairs[pair, 1]] = max(estimate)
      
      
    }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
      Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
      Yj = as.numeric(loci_table[pairs[pair, 2],])
      
      estimate = NULL
      for(haplotype_i in 1:ncol(Yi)){
        estimate = c(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj, freq_table, markers))
      }
      
      relatedness_matrix[pairs[pair, 2],pairs[pair, 1]] = max(estimate)
      
      
    }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
      
      Yi = as.numeric(loci_table[pairs[pair, 1],])
      Yj = as.numeric(loci_table[pairs[pair, 2],])
      
      relatedness_matrix[pairs[pair, 2],pairs[pair, 1]] = fx_get_relatedness(Yi, Yj, freq_table, markers)
      
    }
    
    print(paste("pair", pair, "of", nrow(pairs), sep = " "))
    
  }
  
  return(relatedness_matrix)
  
}
