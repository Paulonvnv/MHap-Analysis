
fx_haplotypes_respect_to_reference = function(ampseq_object,
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
  library(ape)
  library(Biostrings)
  
  # Call 3D7 genome and gff---
  # Call 3D7 reference genome and its corresponding anotation in the gff file
  reference_gff = read.gff(gff_file)
  reference_genome = readDNAStringSet(fasta_file)
  
  # Call set of MHAPs---
  # Set of amplicon markers, This table contains the list mhap, their chromosome location in the 3D7 strain
  markers = ampseq_object$markers
  
  ## Filter Drug resistant markers---
  
  markers_of_interest = markers[grep(paste(gene_names, collapse = "|"),markers$amplicon),]
  
  ### Rename gene common name to gene_ID---
  
  markers_of_interest[['gene_ids']] = NA
  
  for(gene in 1:length(gene_names)){
    markers_of_interest[grepl(gene_names[gene], markers_of_interest[['amplicon']]),][['gene_ids']] = gene_ids[gene]
  }
  
  # Calculates the start and end position of each drugR marker on the CDS of each gene---
  
  ## Start and end position in 3D7 CDSs---
  markers_of_interest[['start_cds']] = NA
  markers_of_interest[['end_cds']] = NA
  markers_of_interest[['strand']] = NA
  markers_of_interest[['ref_length']] = NA
  
  for(amplicon in markers_of_interest$amplicon){ # for each drugR marker
    
    # Gene where the drugR marker is located
    gene = markers_of_interest[markers_of_interest$amplicon == amplicon,][['gene_ids']]
    
    # Filter the genes attributes of the CDS
    temp_gff = reference_gff[grepl(gene, reference_gff$attributes, ignore.case = TRUE) &
                               reference_gff$type == 'CDS',]
    
    # start counting from 0 taking into account the CDSs (exons) where the drugR marker is located
    start_cds = 0
    
    # check how many CDSs does the gene have
    if(nrow(temp_gff) > 1){ # if the gene has 2 or more CDSs then
      
      # Identify on which CDS the drugR marker is located
      cds = which(temp_gff$start < markers_of_interest[markers_of_interest$amplicon== amplicon,]$start &
                    temp_gff$end > markers_of_interest[markers_of_interest$amplicon== amplicon,]$start)
      
      if(cds > 1){ # if the drugR marker is located in CDS 2 or above then
        
        for(pre_cds in 1:(cds - 1)){ # sum the length of all previous CDSs
          start_cds = start_cds + temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1
        }
        
        # sum the position on the located CDS
        start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
        end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
        
        
      }else{ # if the drugR marker is located in CDS 1
        
        # Take the position on the located CDS
        
        start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
        end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
        
      }
    }else{# if the gene has 1 CDS then
      
      # Take the position on the CDS 1
      start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[['start']] + 1
      end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
      
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[['strand']])
      
    }
  }
  
  
  # Generate a reference sequence of each gene---
  # Genes located in the negative strand will not be transformed to their reverse complement yet
  
  # Reference sequences from 3D7---
  ref_seqs = NULL
  
  for(gene in unique(markers_of_interest$gene_ids)){
    
    temp_gff = reference_gff[grepl(gene, reference_gff$attributes, ignore.case = T)&
                               reference_gff$type == 'CDS',]
    
    temp_refseq = NULL
    
    if(nrow(temp_gff)>1){
      for(cds in 1:nrow(temp_gff)){
        temp_refseq = paste0(temp_refseq,
                             subseq(reference_genome[grep(temp_gff[cds,][['seqid']], names(reference_genome))],
                                    start = temp_gff[cds,][['start']],
                                    end = temp_gff[cds,][['end']]))
      }
      
      ref_seqs = c(ref_seqs, as.character(temp_refseq))
      
    }else{
      temp_refseq = subseq(reference_genome[grep(temp_gff[['seqid']], names(reference_genome))],
                           start = temp_gff[['start']],
                           end = temp_gff[['end']])
      
      ref_seqs = c(ref_seqs, as.character(temp_refseq))
    }
    
    markers_of_interest[markers_of_interest$gene_ids == gene,'ref_length'] = nchar(temp_refseq)
    
  }
  
  names(ref_seqs) = unique(markers_of_interest$gene_ids)
  ref_seqs = DNAStringSet(ref_seqs)
  
  
  # Convert cigar format to standard mutation nomenclature PMC1867422---
  
  ## Filter drugR markers---
  moi_loci_abd_table = ampseq_object$loci_abd_table[,grep(paste(gene_names, collapse = "|"), colnames(ampseq_object$loci_abd_table))]
  
  ## Remove read abundace---
  moi_loci_abd_table = gsub(":[0-9]+", "", moi_loci_abd_table)
  
  
  moi_loci_dna_table = moi_loci_abd_table
  moi_loci_aa_table = moi_loci_abd_table
  
  
  for(amplicon in colnames(moi_loci_abd_table)){ # For each amplicon in columns
    for(sample in 1:nrow(moi_loci_abd_table)){ # For each sample in rows
      
      locus = moi_loci_abd_table[sample, amplicon] # Get the genotype in the locus
      
      if(is.na(locus)){ # if the locus is NULL complete the cell with NA
        
        moi_loci_dna_table[sample, amplicon] = NA
        moi_loci_aa_table[sample, amplicon] = NA
        
      }else{
        
        clones = unlist(strsplit(locus, '_')) # get all different clones at that locus
        
        dna_clones = NULL
        aa_clones = NULL
        
        for(clone in clones){ # for each clone get alleles of all different SNPs
          
          alleles = paste0(unlist(strsplit(gsub('[ATCGDI= //.]$', '', clone), '[ATCG]')),
                           unlist(strsplit(gsub('^[0-9]+', '', clone), '[0-9]+')))
          
          dna_alleles = NULL
          aa_alleles = NULL
          
          
          if(clone == "."){ # if allele equals to wild type
            
            codons = "."
            
            }else{
              
              codons = data.frame(alleles = gsub('[0-9]','',alleles), cds_position = sapply(alleles, function(allele){
                
                # position in the mhap
                mhap_position = as.integer(gsub('[ATCGDI=//.]+', '', allele))
                
                # variant or nucleotide found
                mhap_variant = gsub('[0-9]+', '', allele)
                
                # calculate position in the CDS
                
                cds_position = markers_of_interest[markers_of_interest$amplicon == amplicon,'start_cds'] + mhap_position - 1
                
              }))
              
              
              # Identify nucleotide in the reference strain
              codons$ref_variant = sapply(1:nrow(codons), function(x){
                as.character(subseq(ref_seqs[which(names(ref_seqs) ==
                  markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                start = codons[x, 'cds_position'],
                                                end = codons[x,'cds_position']))})
              
              # calculate the aminoacid position
              codons$aa_position = ceiling(codons$cds_position/3)
              codons$first_nucleotide = 3*codons$aa_position - 2
              codons$last_nucleotide = 3*codons$aa_position
              
            }
          
        if(length(codons) > 1){
          
          list_of_codons = unique(codons$aa_position)
          
        }else{
            list_of_codons = codons

        }
        
          
          for(codon in list_of_codons){ # for each allele
            
            if(codon == "."){ # if allele equals to wild type
              
              dna_alleles = c(dna_alleles, 'c.(=)')
              aa_alleles = c(aa_alleles, 'p.(=)')
              
            }else{ # if allele different than the wild type
              
              # get the reference amino acid variant
              ref_aa_variant = ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
                                      as.character(Biostrings::translate(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                                    start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                    end = unique(codons[codons$aa_position == codon,][['last_nucleotide']])))),
                                      as.character(Biostrings::translate(reverseComplement(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                                                      start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                                      end = unique(codons[codons$aa_position == codon,][['last_nucleotide']]))))))
              
              
              codon_positons = unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']])
              
              positions_from_ref = codon_positons[!(unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']]) %in% 
                    codons[codons$aa_position == codon,][['cds_position']])]
              
              positions_from_ref_df = NULL
              
              for(positions_from_ref1 in  positions_from_ref){
                positions_from_ref_df = rbind(positions_from_ref_df,
                      data.frame(alleles = as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                         start = positions_from_ref1,
                                                         end = positions_from_ref1)),
                           cds_position = positions_from_ref1))
              }
              
              
              sample_codon = rbind(codons[codons$aa_position == codon,c('alleles', 'cds_position')],
                                   positions_from_ref_df
              )
              
              sample_codon = sample_codon[order(sample_codon$cds_position),]
              
              mhap_aa_variant = as.character(Biostrings::translate(DNAString(paste0(sample_codon$alleles, collapse = ""))))
        
              
              # 
              # 
              # mhap_aa_variant = if(cds_position - first_nucleotide == 0){
              #   
              #   ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
              #          as.character(translate(DNAString(paste0(
              #            mhap_variant,
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide + 1,
              #                                end = first_nucleotide + 2)))))),
              #          as.character(translate(reverseComplement(DNAString(paste0(
              #            mhap_variant,
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide + 1,
              #                                end = first_nucleotide + 2)))))))
              #   )
              #   
              # }else if(cds_position - first_nucleotide == 1){
              #   
              #   ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
              #          as.character(translate(DNAString(paste0(
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide,
              #                                end = first_nucleotide)),
              #            mhap_variant,
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide + 2,
              #                                end = first_nucleotide + 2))
              #          )))),
              #          as.character(translate(reverseComplement(DNAString(paste0(
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide,
              #                                end = first_nucleotide)),
              #            mhap_variant,
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide + 2,
              #                                end = first_nucleotide + 2))
              #          )))))
              #          
              #   )
              #   
              #   
              # }else if(cds_position - first_nucleotide == 2){
              #   ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
              #          as.character(translate(DNAString(paste0(
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide,
              #                                end = first_nucleotide + 1)),
              #            mhap_variant)))),
              #          as.character(translate(reverseComplement(DNAString(paste0(
              #            as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
              #                                start = first_nucleotide,
              #                                end = first_nucleotide + 1)),
              #            mhap_variant)))))
              #   )
              #   
              # }
            
              
              
              
              
              if(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+"){
                dna_alleles = c(dna_alleles, paste(paste0('c.',codons[codons$aa_position == codon,][['cds_position']],
                                                          codons[codons$aa_position == codon,][['ref_variant']],
                                                          '>',
                                                          codons[codons$aa_position == codon,][['alleles']]), collapse = " "))
                aa_alleles = c(aa_alleles, paste0(ref_aa_variant, codon, mhap_aa_variant))
              }else{
                dna_alleles = c(dna_alleles, paste(paste0('c.',
                                                          markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 
                                                            codons[codons$aa_position == codon,][['cds_position']] + 1,
                                                          codons[codons$aa_position == codon,][['ref_variant']], '>',
                                                          codons[codons$aa_position == codon,][['alleles']]), collapse = ' '))
                aa_alleles = c(aa_alleles, paste0(ref_aa_variant,
                                                  ceiling((markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 3*codon + 1)/3),
                                                  mhap_aa_variant))
                
              }
              
              
            }
            
            
          }
          
          dna_clones[[clone]] = paste(dna_alleles, collapse = " ")
          aa_clones[[clone]] = paste(aa_alleles, collapse = " ")
          
        }
        
        moi_loci_dna_table[sample, amplicon] = paste(dna_clones, collapse = " / ")
        moi_loci_aa_table[sample, amplicon] = paste(aa_clones, collapse = " / ")
        
        
      }
    }
  }
  
  haplotypes_respect_to_reference = list(markers_of_interest = markers_of_interest,
                                         loci_dna_table = moi_loci_dna_table,
                                         loci_aa_table = moi_loci_aa_table,
                                         ref_seqs = ref_seqs)
  
  return(haplotypes_respect_to_reference)
  
}
