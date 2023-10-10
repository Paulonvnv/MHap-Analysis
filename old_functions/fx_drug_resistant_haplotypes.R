

fx_drug_resistant_haplotypes = function(ampseq_object,
                                        reference_alleles = 'drugR_alleles.csv',
                                        gene_names = c('PfDHFR',
                                                       'PfMDR1',
                                                       'PfDHPS',
                                                       'PfKelch13',
                                                       'PF3D7_1447900'),
                                        gene_ids = c('PF3D7_0417200',
                                                     'PF3D7_0523000',
                                                     'PF3D7_0810800',
                                                     'PF3D7_1343700',
                                                     'PF3D7_1447900'),
                                        gff_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff",
                                        fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta",
                                        variables = c('samples', 'Population', 'quarter_of_collection'),
                                        filters = c('Population;Buenaventura,Quibdo,Guapi',
                                                    'quarter_of_collection;2020-Q4,2021-Q1,2021-Q2,2021-Q3,2021-Q4'),
                                        na.var.rm = FALSE){
  
  
  # Call reference alleles
  drugR_reference_alleles = read.csv(reference_alleles)
  
  drugR_reference_alleles['reference'] = gsub('[0-9]+[A-Z]$', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['mutant'] = gsub('^[A-Z][0-9]+', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['position'] = gsub('^[A-Z]', '', gsub('[A-Z]$', '', drugR_reference_alleles$Mutation))
  
  # Define haplotypes respect to a reference genome
  source('fx_haplotypes_respect_to_reference.R')
  
  haplotypes_respect_to_reference = fx_haplotypes_respect_to_reference(ampseq_object,
                                                                       gene_names,
                                                                       gene_ids,
                                                                       gff_file,
                                                                       fasta_file)
  
  
  loci_aa_table = haplotypes_respect_to_reference[['loci_aa_table']]
  
  drug_markers = haplotypes_respect_to_reference$markers_of_interest
  
  ref_seqs = haplotypes_respect_to_reference$ref_seqs
  
  # Empty table to fill cigar outputs
  aacigar_table = matrix(NA,
                         nrow = nrow(loci_aa_table),
                         ncol = length(unique(drugR_reference_alleles$Gene_Id)),
                         dimnames = list(rownames(loci_aa_table),
                                         unique(drugR_reference_alleles$Gene_Id)))
  
  phenotype_table = matrix(NA,
                           nrow = nrow(loci_aa_table),
                           ncol = length(unique(drugR_reference_alleles$Gene_Id)),
                           dimnames = list(rownames(loci_aa_table),
                                           unique(drugR_reference_alleles$Gene_Id)))
  
  for(gene in unique(drugR_reference_alleles$Gene_Id)){ # For each gene
    
    for(amplicon in drug_markers[drug_markers$gene_ids == gene, 'amplicon']){ # for each amplicon in the gene
      
      if(amplicon %in% colnames(loci_aa_table)){
        
        found_positions = loci_aa_table[,amplicon]
        
        found_positions = as.character(unlist(strsplit(as.character(unlist(strsplit(found_positions, " / "))), ' ')))
        
        
        found_positions = unique(as.character(gsub('^[A-Z]',
                                                   '',
                                                   gsub('[A-Z]$',
                                                        '',
                                                        found_positions))))
        
        found_positions = found_positions[!grepl('p',found_positions)]
        found_positions = found_positions[!is.na(found_positions)]
        
        # logic vector of amplified amplicon
        tested_amplicon = drug_markers$amplicon %in% amplicon
        
        for(sample in rownames(loci_aa_table)){ # For each sample
          # Identify all polymorphic positions in reference table
          
          positions = unique(c(found_positions, drugR_reference_alleles[drugR_reference_alleles$Gene_Id == gene, 'position']))
          
          # Identify tested positions
          
          if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
            
            tested_positions = seq(
              floor(
                drug_markers[drug_markers$amplicon == amplicon,
                             'start_cds']/3),
              ceiling(drug_markers[drug_markers$amplicon == amplicon,
                                   'end_cds']/3), 1)
            
          }else{
            
            tested_positions = seq(
              floor(
                (drug_markers[drug_markers$amplicon == amplicon,
                              'ref_length'] - 
                   drug_markers[drug_markers$amplicon == amplicon,
                                'end_cds'] + 1)/3),
              ceiling((drug_markers[drug_markers$amplicon == amplicon,
                                    'ref_length'] - 
                         drug_markers[drug_markers$amplicon == amplicon,
                                      'start_cds'] + 1)/3), 1)
            
          }
          
          
          # Keep polymorphic positions in the reference table that have been tested
          positions = positions[as.integer(positions) %in% tested_positions]
          
          # Select amplicons for defined sample and gene
          sample_clones = loci_aa_table[sample,
                                        amplicon]
          
          
          if(!is.na(sample_clones)){
            
            clones = unlist(strsplit(sample_clones, ' / '))
            
            
            clone_alleles = NULL
            
            for(clone in clones){
              
              if(clone != 'p.(=)'){
                
                sample_positions = gsub('^[A-Z]',
                                        '',
                                        gsub('[A-Z]$',
                                             '',
                                             unlist(strsplit(clone, ' '))))
                
                sample_alleles = gsub('^[A-Z][0-9]+',
                                      '',
                                      unlist(strsplit(clone, ' ')))
                
                remaining_positions = positions[!(positions %in% sample_positions)]
                
                for(position in remaining_positions){
                  
                  sample_positions = c(sample_positions, position)
                  
                  if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           reverseComplement(
                                             subseq(
                                               ref_seqs[grep(gene,names(ref_seqs))],
                                               start = drug_markers[drug_markers$amplicon == amplicon,
                                                                    'ref_length'] - as.numeric(position)*3 + 1,
                                               end = drug_markers[drug_markers$amplicon == amplicon,
                                                                  'ref_length'] - as.numeric(position)*3 + 3)))))
                    
                  }
                }
                
              }else{
                
                sample_positions = NULL
                sample_alleles = NULL
                
                for(position in positions){
                  
                  sample_positions = c(sample_positions, position)
                  
                  if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           reverseComplement(
                                             subseq(
                                               ref_seqs[grep(gene,names(ref_seqs))],
                                               start = drug_markers[drug_markers$amplicon == amplicon,
                                                                    'ref_length'] - as.numeric(position)*3 + 1,
                                               end = drug_markers[drug_markers$amplicon == amplicon,
                                                                  'ref_length'] - as.numeric(position)*3 + 3)))))
                    
                  }
                }
                
              }
              
              clone_alleles = rbind(clone_alleles, sample_alleles)
              
            }
            
            if(length(clones) > 1){
              sample_alleles = apply(clone_alleles, 2, function(x) ifelse(x[1] == x[2], x[1], paste(x[1], x[2], sep = '|')))
              
            }else{
              
              sample_alleles = clone_alleles
              
            }
            
            
          }else{
            
            sample_positions = NULL
            sample_alleles = NULL
            
            for(position in positions){
              
              sample_positions = c(sample_positions, position)
              
              sample_alleles = c(sample_alleles, '?')
              
            }
            
          }
          
          reference_alleles = NULL
          
          for(position in sample_positions){
            
            if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                            drugR_reference_alleles[['position']] == position,]) != 0){
              
              reference_alleles = c(reference_alleles, unique(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                                        drugR_reference_alleles[['position']] == position,'reference']))
            }else{
              
              if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                
                reference_alleles = c(reference_alleles,
                                      as.character(
                                        Biostrings::translate(
                                          subseq(
                                            ref_seqs[grep(gene,names(ref_seqs))],
                                            start = as.numeric(position)*3 - 2,
                                            end = as.numeric(position)*3))))
                
              }else{
                
                reference_alleles = c(reference_alleles,
                                      as.character(
                                        Biostrings::translate(
                                          reverseComplement(
                                            subseq(
                                              ref_seqs[grep(gene,names(ref_seqs))],
                                              start = drug_markers[drug_markers$amplicon == amplicon,
                                                                   'ref_length'] - as.numeric(position)*3 + 1,
                                              end = drug_markers[drug_markers$amplicon == amplicon,
                                                                 'ref_length'] - as.numeric(position)*3 + 3)))))
                
              }
              
              
            }
            
          }
          
          aacigar_haplotype = data.frame(position = as.integer(sample_positions),
                                         aacigar_haplotype = paste(reference_alleles, sample_positions, sample_alleles, sep = ""))
          
          aacigar_haplotype[['aacigar_haplotype']] = ifelse(aacigar_haplotype[['position']] %in% drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene,'position'],
                                                            aacigar_haplotype[['aacigar_haplotype']],
                                                            tolower(aacigar_haplotype[['aacigar_haplotype']]))
          
          sample_phenotype = NULL
          
          for(position in sort(aacigar_haplotype[['position']])){
            
            if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                            drugR_reference_alleles[['position']] == position,]) != 0){
              
              sample_allele = gsub('[A-Z][0-9]+', '', aacigar_haplotype[aacigar_haplotype[['position']] == position, 'aacigar_haplotype'])
              
              reference_sensitive_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                      drugR_reference_alleles[['position']] == position,'reference']
              
              if(sample_allele %in% reference_sensitive_alleles){
                
                sample_phenotype = c(sample_phenotype, 'Sensitive phenotype')  
                
              }else{
                
                reference_resistant_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                        drugR_reference_alleles[['position']] == position,'mutant']
                
                if(sum(strsplit(sample_allele, "\\|")[[1]] %in% reference_resistant_alleles) >= 1){
                  
                  if(grepl('|', sample_allele)){
                    
                    sample_allele = strsplit(sample_allele, "\\|")[[1]]
                    
                    sample_allele = sample_allele[sample_allele %in% reference_resistant_alleles]
                    
                  }
                  
                  phenotype = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                        drugR_reference_alleles[['position']] == position  &
                                                        drugR_reference_alleles[['mutant']] == sample_allele,'Anotation']
                  
                  sample_phenotype = c(sample_phenotype, phenotype)
                  
                }else{
                  
                  if(sample_allele != '?'){
                    
                    sample_phenotype = c(sample_phenotype, paste0(sample_allele, ' variant unreported for position ', position))
                    
                  }else{
                    
                    sample_phenotype = c(sample_phenotype, paste0(amplicon, ' amplicon did not amplify'))
                    
                  }
                  
                }
                
              }
              
            }else{
              
              reference_allele = gsub('[0-9]+[A-Z]', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
              
              alternative_allele = gsub('[A-Z][0-9]+', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
              
              if(reference_allele != alternative_allele){
                
                sample_phenotype = c(sample_phenotype, paste0(aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ' unreported polymorphims respect to 3D7'))  
                
              }
              
            }
            
          }
          
          sample_phenotype = paste(sample_phenotype, collapse = '; ')
          
          aacigar_haplotype = paste(aacigar_haplotype[order(aacigar_haplotype$position),'aacigar_haplotype'], collapse = " ")
          
          
          # Update for each amplicon in the corresponding gene
          if(is.na(aacigar_table[sample, gene])){
            
            aacigar_table[sample, gene] = aacigar_haplotype
            phenotype_table[sample, gene] = sample_phenotype
            
          }else{
            
            aacigar_table[sample, gene] = paste(aacigar_table[sample, gene], aacigar_haplotype, sep = ' ')
            phenotype_table[sample, gene] = paste(phenotype_table[sample, gene], sample_phenotype, sep = '; ')
            
          }
          
        }
        
      }
      
    }
  }
  
  
  mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
  
  poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]
  
  poly1 = gsub('\\|([a-z])', '', poly_aacigar_table, ignore.case = T)
  
  poly2 = gsub('([a-z]\\|)', '', poly_aacigar_table, ignore.case = T)
  
  extended_aacigar_table = rbind(mon_aacigar_table,
                                 poly1,
                                 poly2)
  
  extended_aacigar_table = data.frame(samples = rownames(extended_aacigar_table), extended_aacigar_table)
  
  metadata = ampseq_object$metadata[,variables]
  
  names(metadata) = c('samples', 'var1', 'var2')
  
  extended_aacigar_table = merge(extended_aacigar_table,
                                 metadata,
                                 by = 'samples',
                                 all.x = TRUE)
  
  extended_aacigar_table %<>% pivot_longer(cols = all_of(gene_ids),
                                           names_to = 'gene_ids',
                                           values_to = 'haplotype')
  
  extended_aacigar_table$gene_names = NA
  
  for(gene in 1:length(gene_ids)){
    extended_aacigar_table[extended_aacigar_table$gene_ids == gene_ids[gene],][['gene_names']] = gene_names[gene]
    
  }
  
  
  if(na.var.rm){
    extended_aacigar_table = extended_aacigar_table[!is.na(extended_aacigar_table$var1)|
                                                      (!is.na(extended_aacigar_table$var2)|
                                                         grepl('NA',extended_aacigar_table$var2)),]
  }else{
    extended_aacigar_table %<>% mutate(var1 = case_when(
      is.na(var1) ~ paste(variables[2], 'missing'),
      !is.na(var1) ~ var1),
      var2 = case_when(
        is.na(var2)|grepl('NA',var2) ~ paste(variables[3], 'missing'),
        (!is.na(var2))&(!grepl('NA',var2)) ~ var2))
  }
  
  
  samples_pop_quarter = extended_aacigar_table %>% group_by(var1, var2)%>%
    summarise(count = nlevels(as.factor(samples))) 
  
  haplotype_counts = extended_aacigar_table %>% group_by(gene_names, var1, var2, haplotype)%>%
    summarise(count = n())
  
  if(!is.null(filters)){
    filters = strsplit(filters,';')
    for(temp_filter in 1:length(filters)){
      if(which(variables == filters[[temp_filter]][1]) == 2){
        haplotype_counts %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        samples_pop_quarter %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
      }else if(which(variables == filters[[temp_filter]][1]) == 3){
        haplotype_counts %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        samples_pop_quarter %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
      }
    }
  }
  
  
  haplotype_counts$freq = NA
  
  for(gene in levels(as.factor(haplotype_counts$gene_names))){
    for(Pop in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene,][['var1']]))){
      for(date in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene&
                                                       haplotype_counts$var1 == Pop,][['var2']]))){
        
        haplotype_counts[haplotype_counts$gene_names == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq']] = haplotype_counts[haplotype_counts$gene_names == gene&
                                                                                                              haplotype_counts$var1 == Pop&
                                                                                                              haplotype_counts$var2 == date,][['count']]/
          sum(haplotype_counts[haplotype_counts$gene_names == gene&
                                 haplotype_counts$var1 == Pop&
                                 haplotype_counts$var2 == date,][['count']])
          # samples_pop_quarter[samples_pop_quarter$var1 == Pop&
          #                       samples_pop_quarter$var2 == date,][['count']]
        
      }
    }
  }
  
  haplotype_counts %<>% mutate(gene_haplo = paste(gene_names, haplotype, sep = ": "))
  
  nhaplo = haplotype_counts %>% group_by(gene_names)%>% summarise(nhaplo = nlevels(as.factor(haplotype))) 
  
  max_ncolors = brewer.pal(max(nhaplo$nhaplo), 'Set3')
  
  colors = NULL
  
  for(n in 1:nrow(nhaplo)){
    colors = c(colors, max_ncolors[1:nhaplo[n,][['nhaplo']]])
    
  }
  
  haplo_freq_plot = haplotype_counts%>%
    ggplot(aes(y = freq, x = var2, fill  = gene_haplo)) +
    geom_col()+
    facet_grid(var1 ~ gene_names)+
    scale_fill_manual(values = colors)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    labs(y = 'Frequency in population',
         x = 'Date of Collection',
         fill = 'Gene: Haplotype')
  
  
  names(haplotype_counts) = c(names(haplotype_counts)[1], variables[2:3], names(haplotype_counts)[-1:-3])
  
  
  drug_resistant_haplotypes = list(aacigar_table = aacigar_table,
                                   phenotype_table = phenotype_table,
                                   haplotype_counts =haplotype_counts,
                                   haplo_freq_plot = haplo_freq_plot)
  
  return(drug_resistant_haplotypes)
  
}






















