
drug_loci_aa_table

pf_drug_markers[,c('amplicon', 'Pf3D7_ID')]


gene = "PF3D7_0417200" # DHFR
gene = "PF3D7_0523000" # MDR1
gene = "PF3D7_0810800" # DHPS
gene = "PF3D7_1447900" # MDR2
gene = "PF3D7_1343700" # K13

amplicon = 'PfDHFR_2'

sample = 'SP0126288181'

aacigar_table['SP0101254157',"PF3D7_0417200"]

# Empty table to fill cigar outputs
aacigar_table = matrix(NA,
       nrow = nrow(drug_loci_aa_table),
       ncol = length(unique(drugR_reference_alleles$Gene_Id)),
       dimnames = list(rownames(drug_loci_aa_table),
                       unique(drugR_reference_alleles$Gene_Id)))

phenotype_table = matrix(NA,
                       nrow = nrow(drug_loci_aa_table),
                       ncol = length(unique(drugR_reference_alleles$Gene_Id)),
                       dimnames = list(rownames(drug_loci_aa_table),
                                       unique(drugR_reference_alleles$Gene_Id)))

for(gene in unique(drugR_reference_alleles$Gene_Id)){ # For each gene
  
  for(amplicon in pf_drug_markers[pf_drug_markers$Pf3D7_ID == gene, 'amplicon']){
    
    if(amplicon %in% colnames(drug_loci_aa_table)){
      
      found_positions = drug_loci_aa_table[,
                                           amplicon]
      
      found_positions = as.character(unlist(strsplit(as.character(unlist(strsplit(found_positions, " / "))), ' ')))
      
      
      found_positions = unique(as.character(gsub('^[A-Z]',
                                                 '',
                                                 gsub('[A-Z]$',
                                                      '',
                                                      found_positions))))
      
      found_positions = found_positions[!grepl('p',found_positions)]
      found_positions = found_positions[!is.na(found_positions)]
      
      # logic vector of amplified amplicon
      tested_amplicon = pf_drug_markers$amplicon %in% amplicon
      
      for(sample in rownames(drug_loci_aa_table)){ # For each sample
        # Identify all polymorphic positions in reference table
        
        positions = unique(c(found_positions, drugR_reference_alleles[drugR_reference_alleles$Gene_Id == gene, 'position']))
        
        # Identify tested positions
        
        if(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'][1] == '+'){
          
          tested_positions = seq(
            floor(
              pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                              'start_cds_3D7']/3),
              ceiling(pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                    'end_cds_3D7']/3), 1)
          
        }else{
          
          tested_positions = seq(
            floor(
              (pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                               'length_3D7'] - 
                 pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                 'end_cds_3D7'] + 1)/3),
            ceiling((pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                     'length_3D7'] - 
                       pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                       'start_cds_3D7'] + 1)/3), 1)
          
        }
        
        
        # Keep polymorphic positions in the reference table that have been tested
        positions = positions[as.integer(positions) %in% tested_positions]
        
        # Select amplicons for defined sample and gene
        sample_clones = drug_loci_aa_table[sample,
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
                
                if(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                  
                  sample_alleles = c(sample_alleles,
                                     as.character(
                                       translate(
                                         subseq(
                                           ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                           start = as.numeric(position)*3 - 2,
                                           end = as.numeric(position)*3))))
                  
                }else{
                  
                  sample_alleles = c(sample_alleles,
                                     as.character(
                                       translate(
                                         reverseComplement(
                                           subseq(
                                             ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                             start = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                                     'length_3D7'] - as.numeric(position)*3 + 1,
                                             end = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                                   'length_3D7'] - as.numeric(position)*3 + 3)))))
                  
                }
              }
              
            }else{
              
              sample_positions = NULL
              sample_alleles = NULL
              
              for(position in positions){
                
                sample_positions = c(sample_positions, position)
                
                if(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                  
                  sample_alleles = c(sample_alleles,
                                     as.character(
                                       translate(
                                         subseq(
                                           ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                           start = as.numeric(position)*3 - 2,
                                           end = as.numeric(position)*3))))
                  
                }else{
                  
                  sample_alleles = c(sample_alleles,
                                     as.character(
                                       translate(
                                         reverseComplement(
                                           subseq(
                                             ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                             start = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                                     'length_3D7'] - as.numeric(position)*3 + 1,
                                             end = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                                   'length_3D7'] - as.numeric(position)*3 + 3)))))
                  
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
            
            if(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'][1] == '+'){
              
              reference_alleles = c(reference_alleles,
                                 as.character(
                                   translate(
                                     subseq(
                                       ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                       start = as.numeric(position)*3 - 2,
                                       end = as.numeric(position)*3))))
              
            }else{
              
              reference_alleles = c(reference_alleles,
                                 as.character(
                                   translate(
                                     reverseComplement(
                                       subseq(
                                         ref_seqs_3D7[grep(gene,names(ref_seqs_3D7))],
                                         start = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                                 'length_3D7'] - as.numeric(position)*3 + 1,
                                         end = pf_drug_markers[pf_drug_markers$amplicon == amplicon,
                                                               'length_3D7'] - as.numeric(position)*3 + 3)))))
              
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

write.csv(x = aacigar_table, file = 'aacigar_table.csv', row.names = T, col.names = T)
write.csv(x = phenotype_table, file = 'phenotype_table.csv', row.names = T, col.names = T)



















