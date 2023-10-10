fx_pairwise_relatedness = function(loci_object, parallel = TRUE, w = 1, n = 1){
  library(parallel)
  library(doMC)
  library(svMisc)
  source('estimate_r_and_k.R')
  
  loci_table = loci_object[["loci_table"]]
  freq_table = loci_object[["freq_table"]]
  markers = loci_object[["markers"]]
  
  pairs = as.data.frame(t(combn(rownames(loci_table), 2)))
  
  s = round(seq(1,nrow(pairs)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  pairs = pairs[low:high,]
  
  fx_get_relatedness = function(Yi, Yj, freq_table, markers){
    
    Ys = cbind(Yi, Yj)
    
    rownames(Ys) = colnames(loci_table)
    
    Ys =Ys[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    fs = freq_table[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers = markers[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers[["distance"]] = Inf
    
    for(chromosome in levels(as.factor(no_na_markers[["chromosome"]]))){
      for(amplicon in 1:(nrow(no_na_markers[no_na_markers[["chromosome"]] == chromosome,])-1)){
        no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "distance"] = 
          no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "pos"]
      }
    }
    
    ds = no_na_markers$distance
    
    estimate = estimate_r_and_k(fs = fs,
                                ds = ds,
                                Ys = Ys,
                                warn_fs = F)
    
    return(estimate)
  }
  
  if(parallel){
    registerDoMC(detectCores())
    pairwise_df = foreach(pair = 1:nrow(pairs), .combine = 'rbind') %dopar% {
      
      # If both samples are polyclonal
      if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          for(haplotype_j in 1:ncol(Yj)){
            estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj[,haplotype_j], freq_table, markers))
          }
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the second sample is polyclonal 
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_j in 1:ncol(Yj)){
          estimate = rbind(estimate, fx_get_relatedness(Yi, Yj[,haplotype_j], freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the first sample is polyclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj, freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If both samples are monoclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = fx_get_relatedness(Yi, Yj, freq_table, markers)
        
      }
      
      data.frame(Yi = pairs[pair, 1], Yj = pairs[pair, 2], t(estimate))
      
    }
    
  }else{
    
    pairwise_df = NULL
    
    for(pair in 1:nrow(pairs)){
      
      # If both samples are polyclonal
      if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          for(haplotype_j in 1:ncol(Yj)){
            estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj[,haplotype_j], freq_table, markers))
          }
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the second sample is polyclonal 
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_j in 1:ncol(Yj)){
          estimate = rbind(estimate, fx_get_relatedness(Yi, Yj[,haplotype_j], freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the first sample is polyclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj, freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If both samples are monoclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = fx_get_relatedness(Yi, Yj, freq_table, markers)
        
      }
      
      progress(round(100*pair/nrow(pairs)))
      
      pairwise_df = rbind(pairwise_df, data.frame(Yi = pairs[pair, 1], Yj = pairs[pair, 2], t(estimate)))
      
    }
    
  }
  
  return(relatedness_matrix)
  
}



