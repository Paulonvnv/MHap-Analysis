
library(Biostrings)

# Call 3D7 genome and gff----
# Call 3D7 reference genome and its corresponding anotation in the gff file
Pf3D7_gff = read.gff("reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff")
Pf3D7_genome = readDNAStringSet("reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta")

# Call Dd2 genome and gff----
# Call Dd2 reference genome and its corresponding anotation in the gff file
PfDd2_gff = read.gff("reference/Dd2/PlasmoDB-59_PfalciparumDd2.gff")
PfDd2_genome = readDNAStringSet("reference/Dd2/PlasmoDB-59_PfalciparumDd2_Genome.fasta")

# Call set of MHAPs----
# Set of amplicon markers, This table contains the list mhap, their chromosome location in the 3D7 strain
markers = read.csv("markers.csv")

## Filter Drug resistant markers----
genes = c('PfDHFR', 'PfMDR1', 'PfDHPS', 'PfKelch13', 'PF3D7_1447900')
pf_drug_markers = markers[grep(paste(genes, collapse = "|"),markers$amplicon),]

### Rename gene common name to gene_ID----
pf_drug_markers %<>% mutate(Pf3D7_ID = case_when(
  grepl('DHFR', amplicon) ~ 'PF3D7_0417200',
  grepl('MDR1', amplicon) ~ 'PF3D7_0523000',
  grepl('DHPS', amplicon) ~ 'PF3D7_0810800',
  grepl('Kelch13', amplicon) ~ 'PF3D7_1343700',
  grepl('PF3D7_1447900', amplicon) ~ 'PF3D7_1447900'),
  PfDd2_ID = case_when(
    grepl('DHFR', amplicon) ~ 'PfDd2_040022400',
    grepl('MDR1', amplicon) ~ 'PfDd2_050027900',
    grepl('DHPS', amplicon) ~ 'PfDd2_080015900',
    grepl('Kelch13', amplicon) ~ 'PfDd2_130049500',
    grepl('PF3D7_1447900', amplicon) ~ 'PfDd2_140052300')
  )

# Calculates the start and end position of each drugR marker on the CDS of each gene----

## Start and end position in 3D7 CDSs----
pf_drug_markers[['start_cds_3D7']] = NA
pf_drug_markers[['end_cds_3D7']] = NA
pf_drug_markers[['strand']] = NA
pf_drug_markers[['length_3D7']] = NA
pf_drug_markers[['length_Dd2']] = NA

for(amplicon in pf_drug_markers$amplicon){ # for each drugR marker
  
  # Gene where the drugR marker is located
  gene = pf_drug_markers[pf_drug_markers$amplicon == amplicon,][['Pf3D7_ID']]
  
  # Filter the genes attributes of the CDS
  temp_gff = Pf3D7_gff[grepl(gene, Pf3D7_gff$attributes, ignore.case = T) &
                         Pf3D7_gff$type == 'CDS',]
  
  # start counting from 0 taking into account the CDSs (exons) where the drugR marker is located
  start_cds = 0
  
  # check how many CDSs does the gene have
  if(nrow(temp_gff) > 1){ # if the gene has 2 or more CDSs then
    
    # Identify on which CDS the drugR marker is located
    cds = which(temp_gff$start < pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$start &
            temp_gff$end > pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$start)
    
    if(cds > 1){ # if the drugR marker is located in CDS 2 or above then
      
      for(pre_cds in 1:(cds - 1)){ # sum the length of all previous CDSs
        start_cds = start_cds + temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1
      }
      
      # sum the position on the located CDS
      start_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
      end_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$length - 1
      
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['start_cds_3D7']] = start_cds
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['end_cds_3D7']] = end_cds
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
      
      
    }else{ # if the drugR marker is located in CDS 1
      
      # Take the position on the located CDS
      
      start_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
      end_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$length - 1
      
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['start_cds_3D7']] = start_cds
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['end_cds_3D7']] = end_cds
      pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
      
    }
  }else{# if the gene has 1 CDS then
    
    # Take the position on the CDS 1
    start_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$start - temp_gff[['start']] + 1
    end_cds = start_cds + pf_drug_markers[pf_drug_markers$amplicon== amplicon,]$length - 1
    
    pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['start_cds_3D7']] = start_cds
    pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['end_cds_3D7']] = end_cds
    pf_drug_markers[pf_drug_markers$amplicon== amplicon,][['strand']] = as.character(temp_gff[['strand']])
    
  }
}

## Start and end position in Dd2 CDSs----

# we just have to edit our calculus in 3D7 based on the indels between 3D7 and Dd2

# print(msa(
#   translate(DNAStringSet(c(as.character(ref_seqs_3D7[2]),
#                            as.character(ref_seqs_Dd2[2])
#   )))
# ), show = 'complete')
# 
# print(msa(
#   DNAStringSet(c(as.character(ref_seqs_3D7[2]),
#                            as.character(ref_seqs_Dd2[2])
#   ))
# ), show = 'complete')

# Only PfMDR1_3 and PfMDR1_4 are affected by indels

pf_drug_markers[['start_cds_Dd2']] = ifelse(pf_drug_markers$amplicon %in% c('PfMDR1_3','PfMDR1_4'), pf_drug_markers[['start_cds_3D7']] + 3, pf_drug_markers[['start_cds_3D7']])
pf_drug_markers[['end_cds_Dd2']] = ifelse(pf_drug_markers$amplicon %in% c('PfMDR1_3','PfMDR1_4'), pf_drug_markers[['end_cds_3D7']] + 3, pf_drug_markers[['end_cds_3D7']])


# Generate a reference sequence of each gene----
# Genes located in the negative strand will not be transformed to their reverse complement yet

# Reference sequences from 3D7----
ref_seqs_3D7 = NULL

for(gene in unique(pf_drug_markers$Pf3D7_ID)){
  
  temp_gff = Pf3D7_gff[grepl(gene, Pf3D7_gff$attributes, ignore.case = T)&
                         Pf3D7_gff$type == 'CDS',]
  
  temp_refseq = NULL
  
  if(nrow(temp_gff)>1){
    for(cds in 1:nrow(temp_gff)){
      temp_refseq = paste0(temp_refseq,
                           subseq(Pf3D7_genome[grep(temp_gff[cds,][['seqid']], names(Pf3D7_genome))],
                                  start = temp_gff[cds,][['start']],
                                  end = temp_gff[cds,][['end']]))
    }
    #if(temp_gff[1,][['strand']] == '+'){
      ref_seqs_3D7 = c(ref_seqs_3D7, as.character(temp_refseq))
    #}else{
     # ref_seqs_3D7 = c(ref_seqs_3D7, as.character(reverseComplement(DNAString(temp_refseq))))
    #}
  }else{
    temp_refseq = subseq(Pf3D7_genome[grep(temp_gff[['seqid']], names(Pf3D7_genome))],
                         start = temp_gff[['start']],
                         end = temp_gff[['end']])
    
    #if(temp_gff[['strand']] == '+'){
      ref_seqs_3D7 = c(ref_seqs_3D7, as.character(temp_refseq))
    #}else{
     # ref_seqs_3D7 = c(ref_seqs_3D7, as.character(reverseComplement(temp_refseq)))
    #}
  }
  
  pf_drug_markers[pf_drug_markers$Pf3D7_ID == gene,'length_3D7'] = nchar(temp_refseq)
  
}

names(ref_seqs_3D7) = unique(pf_drug_markers$Pf3D7_ID)
ref_seqs_3D7 = DNAStringSet(ref_seqs_3D7)

# Reference sequences from Dd2----
ref_seqs_Dd2 = NULL

for(gene in unique(pf_drug_markers$PfDd2_ID)){
  
  temp_gff = PfDd2_gff[grepl(gene, PfDd2_gff$attributes, ignore.case = T)&
                         PfDd2_gff$type == 'CDS',]
  
  temp_refseq = NULL
  
  if(nrow(temp_gff)>1){
    for(cds in 1:nrow(temp_gff)){
      temp_refseq = paste0(temp_refseq,
                           subseq(PfDd2_genome[grep(temp_gff[cds,][['seqid']], names(PfDd2_genome))],
                                  start = temp_gff[cds,][['start']],
                                  end = temp_gff[cds,][['end']]))
    }
    #if(temp_gff[1,][['strand']] == '+'){
    ref_seqs_Dd2 = c(ref_seqs_Dd2, as.character(temp_refseq))
    #}else{
    # ref_seqs_Dd2 = c(ref_seqs_Dd2, as.character(reverseComplement(DNAString(temp_refseq))))
    #}
  }else{
    temp_refseq = subseq(PfDd2_genome[grep(temp_gff[['seqid']], names(PfDd2_genome))],
                         start = temp_gff[['start']],
                         end = temp_gff[['end']])
    
    #if(temp_gff[['strand']] == '+'){
    ref_seqs_Dd2 = c(ref_seqs_Dd2, as.character(temp_refseq))
    #}else{
    # ref_seqs_Dd2 = c(ref_seqs_Dd2, as.character(reverseComplement(temp_refseq)))
    #  }
  }
  
  pf_drug_markers[pf_drug_markers$PfDd2_ID == gene,'length_Dd2'] = nchar(temp_refseq)
  
}

names(ref_seqs_Dd2) = unique(pf_drug_markers$PfDd2_ID)
ref_seqs_Dd2 = DNAStringSet(ref_seqs_Dd2)

# Convert cigar format to standard mutation nomenclature PMC1867422----

## Filter drugR markers----
drug_loci_abd_table = ampseq$loci_abd_table[,grep(paste(genes, collapse = "|"), colnames(ampseq$loci_abd_table))]

## Remove read abundace----
drug_loci_abd_table = gsub(":[0-9]+", "", drug_loci_abd_table)


drug_loci_dna_table = drug_loci_abd_table
drug_loci_aa_table = drug_loci_abd_table


for(amplicon in colnames(drug_loci_abd_table)){ # For each amplicon in columns
  for(sample in 1:nrow(drug_loci_abd_table)){ # For each sample in rows
    
    locus = drug_loci_abd_table[sample, amplicon] # Get the genotype in the locus
    
    if(is.na(locus)){ # if the locus is NULL complete the cell with NA
      
      drug_loci_dna_table[sample, amplicon] = NA
      drug_loci_aa_table[sample, amplicon] = NA
      
    }else{
      
      clones = unlist(strsplit(locus, '_')) # get all different clones at that locus
      
      dna_clones = NULL
      aa_clones = NULL
      
      for(clone in clones){ # for each clone get alleles of all different SNPs
        
        alleles = paste0(unlist(strsplit(gsub('[ATCGDI= //.]$', '', clone), '[ATCG]')),
                         unlist(strsplit(gsub('^[0-9]+', '', clone), '[0-9]+')))
        
        dna_alleles = NULL
        aa_alleles = NULL
        
        for(allele in alleles){ # for each allele
          
          if(allele == "."){ # if allele equals to wild type
            
            dna_alleles = c(dna_alleles, 'c.(=)')
            aa_alleles = c(aa_alleles, 'p.(=)')
            
          }else{ # if allele different than the wild type
            
            # position in the mhap
            mhap_position = as.integer(gsub('[ATCGDI=//.]+', '', allele))
            
            # variant or nucleotide found
            mhap_variant = gsub('[0-9]+', '', allele)
            
            # calculate position in the CDS
            
            cds_position = pf_drug_markers[pf_drug_markers$amplicon == amplicon,'start_cds_3D7'] + mhap_position - 1
            
            # Identify nucleotide in the reference strain
            ref_variant = as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7) ==
                                                                   pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                              start = cds_position,
                                              end = cds_position))
            
            # calculate the aminoacid position
            aa_position = ceiling(cds_position/3)
            
            # calculate the position of the first and last nucleotide of the codon
            first_nucleotide = 3*aa_position - 2
            last_nucleotide = 3*aa_position
            
            # get the reference amino acid variant
            ref_aa_variant = ifelse(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'] == "+",
                                    as.character(translate(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                                                  start = first_nucleotide,
                                                                  end = last_nucleotide))),
                                    as.character(translate(reverseComplement(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                                                                    start = first_nucleotide,
                                                                                    end = last_nucleotide)))))
            
            mhap_aa_variant = if(cds_position - first_nucleotide == 0){
              
              ifelse(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'] == "+",
                     as.character(translate(DNAString(paste0(
                       mhap_variant,
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide + 1,
                                           end = first_nucleotide + 2)))))),
                     as.character(translate(reverseComplement(DNAString(paste0(
                       mhap_variant,
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide + 1,
                                           end = first_nucleotide + 2)))))))
              )
              
            }else if(cds_position - first_nucleotide == 1){
              
              ifelse(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'] == "+",
                     as.character(translate(DNAString(paste0(
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide,
                                           end = first_nucleotide)),
                       mhap_variant,
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide + 2,
                                           end = first_nucleotide + 2))
                     )))),
                     as.character(translate(reverseComplement(DNAString(paste0(
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide,
                                           end = first_nucleotide)),
                       mhap_variant,
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide + 2,
                                           end = first_nucleotide + 2))
                     )))))
                     
              )
              
              
            }else if(cds_position - first_nucleotide == 2){
              ifelse(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'] == "+",
                     as.character(translate(DNAString(paste0(
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide,
                                           end = first_nucleotide + 1)),
                       mhap_variant)))),
                     as.character(translate(reverseComplement(DNAString(paste0(
                       as.character(subseq(ref_seqs_3D7[which(names(ref_seqs_3D7)==pf_drug_markers[pf_drug_markers$amplicon == amplicon,'Pf3D7_ID'])],
                                           start = first_nucleotide,
                                           end = first_nucleotide + 1)),
                       mhap_variant)))))
              )
              
            }
            
            if(pf_drug_markers[pf_drug_markers$amplicon == amplicon,'strand'] == "+"){
              dna_alleles = c(dna_alleles, paste0('c.', cds_position, ref_variant, '>', mhap_variant))
              aa_alleles = c(aa_alleles, paste0(ref_aa_variant, aa_position, mhap_aa_variant))
            }else{
              dna_alleles = c(dna_alleles, paste0('c.',
                                                  pf_drug_markers[pf_drug_markers$amplicon == amplicon,'length_3D7'] - cds_position + 1,
                                                  ref_variant, '>', mhap_variant))
              aa_alleles = c(aa_alleles, paste0(ref_aa_variant,
                                                ceiling((pf_drug_markers[pf_drug_markers$amplicon == amplicon,'length_3D7'] - cds_position + 1)/3),
                                                mhap_aa_variant))
              
            }
            
            
          }
          
          
        }
        
        dna_clones[[clone]] = paste(dna_alleles, collapse = " ")
        aa_clones[[clone]] = paste(aa_alleles, collapse = " ")
        
      }
      
      drug_loci_dna_table[sample, amplicon] = paste(dna_clones, collapse = " / ")
      drug_loci_aa_table[sample, amplicon] = paste(aa_clones, collapse = " / ")
      
      
    }
  }
}


View(drug_loci_aa_table)

# 
# write.csv(drug_loci_aa_table, 'drugR_3D7.csv', row.names = T)
# 
# write.csv(drug_loci_aa_table_dd2, 'drugR_Dd2.csv', row.names = T)
# 


drugR_reference_alleles = read.csv('drugR_alleles.csv')

drugR_reference_alleles['reference'] = gsub('[0-9]+[A-Z]$', '', drugR_reference_alleles$Mutation)
drugR_reference_alleles['mutant'] = gsub('^[A-Z][0-9]+', '', drugR_reference_alleles$Mutation)
drugR_reference_alleles['position'] = gsub('^[A-Z]', '', gsub('[A-Z]$', '', drugR_reference_alleles$Mutation))


#drug_loci_aa_table = gsub('^[A-Z]|\\b[A-Z]', '', drug_loci_aa_table)

