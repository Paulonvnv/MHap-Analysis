
# write_loci----

write_loci = function(loci_object, format = c('excel', 'csv', 'tsv'), name = 'wb.xlsx', sep = '\t'){
  
  if(format == 'excel'){
    
    if(file.exists(name)){
      file.remove(name)
    }
    
    excel_wb = loadWorkbook(name, create = T)
    
    for(temp_slot in c('loci_table', 
                       'metadata', 
                       'markers', 
                       'freq_table',
                       'ampseq_alleles_per_locus')){
      
      if(temp_slot == 'loci_table'){
        print("Printing loci_table slot...")
        
        if(!is.null(slot(loci_object, temp_slot))){
          temp_sheet = data.frame(Sample_id = rownames(slot(loci_object, temp_slot)),
                                  as.data.frame(slot(loci_object, temp_slot)))
        } else {
          temp_sheet = NULL
        } 
      }else if(temp_slot == 'markers'){
        print("Printing markers slot...")
        if(!is.null(slot(loci_object, temp_slot))) {
          temp_sheet = as.data.frame(slot(loci_object, temp_slot))
          
          if(sum(is.infinite(temp_sheet[['distance']])) > 0){
            temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
          }
        } else {
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'freq_table'){
        
        
        temp_sheet = NULL
        
        if(!is.null(slot(loci_object, temp_slot))){
          
          print("Printing freq_table slots...")
          
          temp_freq_table = slot(loci_object, temp_slot)
          
          # write loci_table
          
          if(is.list(temp_freq_table)){
            
            for(temp_sub_slot in names(temp_freq_table)){
              
              writeWorksheet(excel_wb,
                             temp_freq_table[[temp_sub_slot]],
                             sheet = paste0('freq_table_', temp_sub_slot),
                             header = T)
              
            }
            
          }else if(is.data.frame(temp_freq_table)){
            
            writeWorksheet(excel_wb,
                           temp_freq_table,
                           sheet = temp_slot,
                           header = T)
            
          }
          
        }
        
        
      }else if(temp_slot == 'ampseq_alleles_per_locus'){
        
        temp_sheet = NULL
        
        # if(!is.null(slot(loci_object, temp_slot))){
        #   
        #   print("Printing ampseq_alleles_per_locus slots...")
        #   
        #   temp_ampseq_alleles_per_locus = slot(loci_object, temp_slot)
        #   
        #   # write loci_table
        #   
        #   if(is.list(temp_freq_table)){
        #     
        #     for(temp_sub_slot in names(temp_freq_table)){
        #       
        #       write.table(temp_freq_table[[temp_sub_slot]],
        #                   file.path(name, paste0('freq_table', temp_sub_slot, '.tsv')), 
        #                   quote = F, row.names = F, sep = sep)
        #       
        #     }
        #     
        #   }else if(is.data.frame(temp_freq_table)){
        #     
        #     write.csv(temp_freq_table, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F, sep = sep)
        #     
        #   }
        #   
        # }
        
        
      }else{
        
        if(!is.null(slot(loci_object, temp_slot))){
          temp_sheet = as.data.frame(slot(loci_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        createSheet(excel_wb, name = temp_slot)
        
        writeWorksheet(excel_wb,
                       temp_sheet,
                       sheet = temp_slot,
                       header = T)
      }
      
      
    }
    
    saveWorkbook(excel_wb)
    
  }else if(format == 'csv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    dir.create(name)
    
    for(temp_slot in c('loci_table', 
                       'metadata', 
                       'markers', 
                       'freq_table',
                       'ampseq_alleles_per_locus')){
      
      if(temp_slot == 'loci_table'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(loci_object, temp_slot)),
                                as.data.frame(slot(loci_object, temp_slot)))
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(loci_object, temp_slot))
        
        if(sum(is.infinite(temp_sheet[['distance']])) > 0){
          temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        }
        
        
      }else if(temp_slot == 'freq_table'){
        
        
        temp_sheet = NULL
        
        if(!is.null(slot(loci_object, temp_slot))){
          
          print("Printing freq_table slots...")
          
          temp_freq_table = slot(loci_object, temp_slot)
          
          # write loci_table
          
          if(is.list(temp_freq_table)){
            
            for(temp_sub_slot in names(temp_freq_table)){
              
              write.csv(temp_freq_table[[temp_sub_slot]],
                        file.path(name, paste0('freq_table_', temp_sub_slot, '.csv')), 
                        quote = F, row.names = F)
              
            }
            
          }else if(is.data.frame(temp_freq_table)){
            
            write.csv(temp_freq_table, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F)
            
          }
          
        }
        
        
      }else if(temp_slot == 'ampseq_alleles_per_locus'){
        
        temp_sheet = NULL
        
        # if(!is.null(slot(loci_object, temp_slot))){
        #   
        #   print("Printing ampseq_alleles_per_locus slots...")
        #   
        #   temp_ampseq_alleles_per_locus = slot(loci_object, temp_slot)
        #   
        #   # write loci_table
        #   
        #   if(is.list(temp_freq_table)){
        #     
        #     for(temp_sub_slot in names(temp_freq_table)){
        #       
        #       write.table(temp_freq_table[[temp_sub_slot]],
        #                   file.path(name, paste0('freq_table', temp_sub_slot, '.tsv')), 
        #                   quote = F, row.names = F, sep = sep)
        #       
        #     }
        #     
        #   }else if(is.data.frame(temp_freq_table)){
        #     
        #     write.csv(temp_freq_table, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F, sep = sep)
        #     
        #   }
        #   
        # }
        
        
      }else{
        
        if(!is.null(slot(loci_object, temp_slot))){
          temp_sheet = as.data.frame(slot(loci_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.csv(temp_sheet, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F)
        
      }
      
      
    }
    
    
  }else if(format == 'tsv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    dir.create(name)
    
    for(temp_slot in c('loci_table', 
                       'metadata', 
                       'markers', 
                       'freq_table',
                       'ampseq_alleles_per_locus')){
      
      if(temp_slot == 'loci_table'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(loci_object, temp_slot)),
                                as.data.frame(slot(loci_object, temp_slot)))
        
      }else if(temp_slot %in% c('asv_seqs', 'asv_seqs_masked')){
        
        if(!is.null(slot(loci_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(loci_object, temp_slot)),
                                  asv_seq = as.character(slot(loci_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(loci_object, temp_slot))
        
        if(sum(is.infinite(temp_sheet[['distance']])) > 0){
          temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        }
        
        
      }else if(temp_slot == 'freq_table'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(loci_object, temp_slot))){
          
          print("Printing freq_table slots...")
          
          temp_freq_table = slot(loci_object, temp_slot)
          
          # write loci_table
          
          if(is.list(temp_freq_table)){
            
            for(temp_sub_slot in names(temp_freq_table)){
              
              write.table(temp_freq_table[[temp_sub_slot]],
                        file.path(name, paste0('freq_table_', temp_sub_slot, '.tsv')), 
                        quote = F, row.names = F, sep = sep)
              
            }
            
          }else if(is.matrix(temp_freq_table)){
            
            write.table(temp_freq_table, paste0(file.path(name, temp_slot), '.tsv'), quote = F, row.names = T, sep = sep)
            
          }
          
        }
        
        
      }else if(temp_slot == 'ampseq_alleles_per_locus'){
        
        temp_sheet = NULL
        
        # if(!is.null(slot(loci_object, temp_slot))){
        #   
        #   print("Printing ampseq_alleles_per_locus slots...")
        #   
        #   temp_ampseq_alleles_per_locus = slot(loci_object, temp_slot)
        #   
        #   # write loci_table
        #   
        #   if(is.list(temp_freq_table)){
        #     
        #     for(temp_sub_slot in names(temp_freq_table)){
        #       
        #       write.table(temp_freq_table[[temp_sub_slot]],
        #                   file.path(name, paste0('freq_table', temp_sub_slot, '.tsv')), 
        #                   quote = F, row.names = F, sep = sep)
        #       
        #     }
        #     
        #   }else if(is.data.frame(temp_freq_table)){
        #     
        #     write.csv(temp_freq_table, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F, sep = sep)
        #     
        #   }
        #   
        # }
        
        
      }else{
        
        if(!is.null(slot(loci_object, temp_slot))){
          temp_sheet = as.data.frame(slot(loci_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.table(temp_sheet, paste0(file.path(name, temp_slot), '.tsv'), quote = F, row.names = F, sep = sep)
        
      }
      
      
    }
    
    
  }
  
}














# read_loci----

read_loci = function(file = NULL, format = 'tsv', sep = '\t'){
  
  loci_object = create_loci()
  
  if(format == 'excel'){
    
    # temp_wb = loadWorkbook(file)
    # 
    # for(sheet in getSheets(temp_wb)){
    #   print(paste0("Loading sheet ", sheet))
    #   if(sheet == 'loci_table'){
    #     
    #     temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #     temp_sheet_rownames = temp_sheet[,1]
    #     temp_sheet = as.matrix(temp_sheet[,-1])
    #     rownames(temp_sheet) = temp_sheet_rownames
    #     
    #     slot(loci_object, sheet, check = TRUE) = temp_sheet
    #     
    #   }else if(sheet %in% c('metadata')){
    #     
    #     temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #     temp_sheet_rownames = temp_sheet[,1]
    #     rownames(temp_sheet) = temp_sheet_rownames
    #     
    #     slot(loci_object, sheet, check = TRUE) = temp_sheet
    #     
    #   }else if(sheet == 'markers'){
    #     
    #     temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #     if (any(is.na(temp_sheet[['distance']]))){
    #       temp_sheet[['distance']][is.na(temp_sheet[['distance']])] = Inf 
    #     }
    #     slot(loci_object, sheet, check = TRUE) = temp_sheet
    #     
    #   }else if(grepl('freq_table', sheet)){
    #     
    #     if(sheet == 'discarded_loci_gt'){
    #       temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #       temp_sheet_rownames = temp_sheet[,1]
    #       temp_sheet = as.matrix(temp_sheet[,-1])
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       
    #       loci_object@discarded_loci[['gt']] = temp_sheet
    #       
    #     }else if(sheet == 'discarded_loci_markers'){
    #       
    #       temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #       if (any(is.na(temp_sheet[['distance']]))){
    #         temp_sheet[['distance']][is.na(temp_sheet[['distance']])] = Inf 
    #       }
    #       loci_object@discarded_loci[['markers']] = temp_sheet
    #       
    #     }else if(sheet == 'discarded_loci_loci_performance'){
    #       
    #       temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #       temp_sheet_rownames = temp_sheet[,1]
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       loci_object@discarded_loci[['loci_performance']] = temp_sheet
    #       
    #     }
    #     
    #   }else{
    #     
    #     temp_sheet = readWorksheet(temp_wb, sheet = sheet)
    #     slot(loci_object, sheet, check = TRUE) = temp_sheet
    #     
    #   }
    # }
    
  }else if(format == 'csv'){
    
    # for(sheet in list.files(file)){
    #   if(sheet == 'loci_table.csv'){
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     temp_sheet_rownames = temp_sheet[,1]
    #     temp_sheet = as.matrix(temp_sheet[,-1])
    #     rownames(temp_sheet) = temp_sheet_rownames
    #     
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }else if(sheet %in% c('asv_seqs.csv', 'asv_seqs_masked.csv')){
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     temp_sheet_names = temp_sheet[[1]]
    #     temp_sheet = DNAStringSet(temp_sheet[[2]])
    #     names(temp_sheet) = temp_sheet_names
    #     
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }else if(sheet %in% c('metadata.csv', 'loci_performance.csv')){
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     temp_sheet_rownames = temp_sheet[,1]
    #     rownames(temp_sheet) = temp_sheet_rownames
    #     
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }else if(sheet %in% c('vcf_like.csv')){
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     temp_sheet_rownames = paste0(temp_sheet[,1], '_', temp_sheet[,2])
    #     rownames(temp_sheet) = temp_sheet_rownames
    #     
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }else if(sheet == 'markers.csv'){
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     if(sum(is.na(temp_sheet[['distance']]))>0){
    #       temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf  
    #     }
    #     
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }else if(grepl('discarded_loci', sheet)){
    #     
    #     if(sheet == 'discarded_loci_gt.csv'){
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       temp_sheet = as.matrix(temp_sheet[,-1])
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       
    #       loci_object@discarded_loci[['gt']] = temp_sheet
    #       
    #     }else if(sheet == 'discarded_loci_markers.csv'){
    #       
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       if(sum(is.na(temp_sheet[['distance']]))>0){
    #         temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
    #       }
    #       
    #       loci_object@discarded_loci[['markers']] = temp_sheet
    #       
    #     }else if(sheet == 'discarded_loci_loci_performance.csv'){
    #       
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       loci_object@discarded_loci[['loci_performance']] = temp_sheet
    #       
    #     }
    #     
    #   }else if(grepl('discarded_samples', sheet)){
    #     
    #     if(sheet == 'discarded_samples_gt.csv'){
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       temp_sheet = as.matrix(temp_sheet[,-1])
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       
    #       loci_object@discarded_samples[['gt']] = temp_sheet
    #       
    #     }else if(sheet == 'discarded_samples_metadata.csv'){
    #       
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       loci_object@discarded_samples[['metadata']] = temp_sheet
    #       
    #     }
    #     
    #   }else if(grepl('controls', sheet)){
    #     
    #     if(sheet == 'controls_gt.csv'){
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       temp_sheet = as.matrix(temp_sheet[,-1])
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       
    #       loci_object@controls[['gt']] = temp_sheet
    #       
    #     }else if(sheet == 'controls_metadata.csv'){
    #       
    #       temp_sheet = read.csv(file.path(file, sheet))
    #       temp_sheet_rownames = temp_sheet[,1]
    #       rownames(temp_sheet) = temp_sheet_rownames
    #       loci_object@controls[['metadata']] = temp_sheet
    #       
    #     }
    #     
    #   }else{
    #     
    #     temp_sheet = read.csv(file.path(file, sheet))
    #     slot(loci_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
    #     
    #   }
    # }
    
  }else if(format == 'tsv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'loci_table.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), sep = sep, header = T)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(loci_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata.tsv', 'loci_performance.tsv')){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(loci_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        
        if(sum(is.na(temp_sheet[['distance']])) > 0){
          temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
        }
        slot(loci_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(grepl('freq_table', sheet)){
        
        if(sheet == 'freq_table.tsv'){
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep, row.names = 1)
          temp_sheet = as.matrix(temp_sheet)
          
          loci_object@freq_table = temp_sheet
          
        }else{
          
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep, row.names = 1)
          # temp_sheet_rownames = temp_sheet[,1]
          # rownames(temp_sheet) = temp_sheet_rownames
          loci_object@freq_table[[gsub('(freq_table_|.tsv)', '', sheet)]] = temp_sheet
          
        }
        
      }else{
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        slot(loci_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
    if(sum(row.names(loci_object@loci_table) != loci_object@metadata$Sample_id) > 0){
      
      row.names(loci_object@loci_table) = loci_object@metadata$Sample_id
      
    }
    
    
  }else if(format == 'json'){
    # In development
  }
  
  return(loci_object)
  
}

