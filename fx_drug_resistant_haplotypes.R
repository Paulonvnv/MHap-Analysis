

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
                                        fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta"){
  
  drugR_reference_alleles = read.csv(reference_alleles)
  
  drugR_reference_alleles['reference'] = gsub('[0-9]+[A-Z]$', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['mutant'] = gsub('^[A-Z][0-9]+', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['position'] = gsub('^[A-Z]', '', gsub('[A-Z]$', '', drugR_reference_alleles$Mutation))
  
  
  source('fx_haplotypes_respect_to_reference.R')
  
  haplotypes_respect_to_reference = fx_haplotypes_respect_to_reference(ampseq,
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
    
    for(amplicon in drug_markers[drug_markers$gene_ids == gene, 'amplicon']){
      
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
                                         translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         translate(
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
                                         translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         translate(
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
                                        translate(
                                          subseq(
                                            ref_seqs[grep(gene,names(ref_seqs))],
                                            start = as.numeric(position)*3 - 2,
                                            end = as.numeric(position)*3))))
                
              }else{
                
                reference_alleles = c(reference_alleles,
                                      as.character(
                                        translate(
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
  
  drug_resistant_haplotypes = list(aacigar_table = aacigar_table,
                                   phenotype_table = phenotype_table)
  
  return(drug_resistant_haplotypes)
  
}






















