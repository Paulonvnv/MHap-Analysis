full_metadata = read.csv('../AMPLseq_Gates/Colombian_reports/Gates_Colombia_metadata.csv')

run1 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_020_021/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run2 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_035_037/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run3 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_036_038_REDO/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run4 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_039_040_REDO2/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run5 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_045_046/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run6 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_047_048/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run7 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_049_050/dada2/run_dada2/CIGARVariants.out.tsv', header = T)
run8 = read.table('../AMPLseq_Gates/Colombian_reports/colombia/MiSeq_AMPLseq_Gates_051/dada2/run_dada2/CIGARVariants.out.tsv', header = T)


run1$Sample_id = gsub('_S\\d+$', '', run1$sample)
run2$Sample_id = gsub('_S\\d+$', '', run2$sample)
run3$Sample_id = gsub('_S\\d+$', '', run3$sample)
run4$Sample_id = gsub('_S\\d+$', '', run4$sample)
run5$Sample_id = gsub('_S\\d+$', '', run5$sample)
run6$Sample_id = gsub('_S\\d+$', '', run6$sample)
run7$Sample_id = gsub('_S\\d+$', '', run7$sample)
run8$Sample_id = gsub('_S\\d+$', '', run8$sample)

run3[run3$Sample_id  == "SP001254204",][['Sample_id']] = "SP0101254204"


genotyped_samples = c(run1$Sample_id,
                      run2$Sample_id,
                      run3$Sample_id,
                      run4$Sample_id,
                      run5$Sample_id,
                      run6$Sample_id,
                      run7$Sample_id,
                      run8$Sample_id)

genotyped_samples = genotyped_samples[grepl('SP', ignore.case = T, genotyped_samples)]


genotyped_samples = unique(genotyped_samples)


metadata = full_metadata[full_metadata$Sample_id %in% genotyped_samples,]

metadata$Sample_id2 = metadata$Sample_id

metadata$Sample_id = c(paste0('ID0000', 1:9),
                       paste0('ID000', 10:99),
                       paste0('ID00', 100:999),
                       paste0('ID0', 1000:1300))


run1$recoded_Sample_id = NA
for(sample in run1$Sample_id){
  if(sample %in% genotyped_samples){
    run1[run1$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run1[run1$Sample_id == sample,][['recoded_Sample_id']] = run1[run1$Sample_id == sample,][['Sample_id']] 
  }
}


run2$recoded_Sample_id = NA
for(sample in run2$Sample_id){
  if(sample %in% genotyped_samples){
    run2[run2$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run2[run2$Sample_id == sample,][['recoded_Sample_id']] = run2[run2$Sample_id == sample,][['Sample_id']] 
  }
}


run3$recoded_Sample_id = NA
for(sample in run3$Sample_id){
  if(sample %in% genotyped_samples){
    run3[run3$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run3[run3$Sample_id == sample,][['recoded_Sample_id']] = run3[run3$Sample_id == sample,][['Sample_id']] 
  }
}

run4$recoded_Sample_id = NA
for(sample in run4$Sample_id){
  if(sample %in% genotyped_samples){
    run4[run4$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run4[run4$Sample_id == sample,][['recoded_Sample_id']] = run4[run4$Sample_id == sample,][['Sample_id']] 
  }
}

run5$recoded_Sample_id = NA
for(sample in run5$Sample_id){
  if(sample %in% genotyped_samples){
    run5[run5$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run5[run5$Sample_id == sample,][['recoded_Sample_id']] = run5[run5$Sample_id == sample,][['Sample_id']] 
  }
}

run6$recoded_Sample_id = NA
for(sample in run6$Sample_id){
  if(sample %in% genotyped_samples){
    run6[run6$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run6[run6$Sample_id == sample,][['recoded_Sample_id']] = run6[run6$Sample_id == sample,][['Sample_id']] 
  }
}


run7$recoded_Sample_id = NA
for(sample in run7$Sample_id){
  if(sample %in% genotyped_samples){
    run7[run7$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run7[run7$Sample_id == sample,][['recoded_Sample_id']] = run7[run7$Sample_id == sample,][['Sample_id']] 
  }
}

run8$recoded_Sample_id = NA
for(sample in run8$Sample_id){
  if(sample %in% genotyped_samples){
    run8[run8$Sample_id == sample,][['recoded_Sample_id']] = metadata[metadata$Sample_id2  == sample,][['Sample_id']]
  }else{
    run8[run8$Sample_id == sample,][['recoded_Sample_id']] = run8[run8$Sample_id == sample,][['Sample_id']] 
  }
}


run1$recoded_Sample_id
run1 = run1[run1$recoded_Sample_id !="A017",]
run1$sample = paste(run1$recoded_Sample_id, gsub('^.+_','',run1$sample), sep = "_")
run1$Sample_id = NULL
run1$recoded_Sample_id = NULL


run2$recoded_Sample_id
run2$sample = paste(run2$recoded_Sample_id, gsub('^.+_','',run2$sample), sep = "_")
run2$Sample_id = NULL
run2$recoded_Sample_id = NULL


run3$recoded_Sample_id
run3$sample = paste(run3$recoded_Sample_id, gsub('^.+_','',run3$sample), sep = "_")
run3$Sample_id = NULL
run3$recoded_Sample_id = NULL


run4$recoded_Sample_id
run4$sample = paste(run4$recoded_Sample_id, gsub('^.+_','',run4$sample), sep = "_")
run4$Sample_id = NULL
run4$recoded_Sample_id = NULL


run5$recoded_Sample_id
run5[grepl('Control',run5$Sample_id),]$recoded_Sample_id = c('Dd2_pos1', 'Dd2_pos2')
run5$sample = paste(run5$recoded_Sample_id, gsub('^.+_','',run5$sample), sep = "_")
run5$Sample_id = NULL
run5$recoded_Sample_id = NULL


run6$recoded_Sample_id
run6[grepl('Control',run6$Sample_id),]$recoded_Sample_id = c('Dd2_pos1', 'Dd2_pos2')
run6$sample = paste(run6$recoded_Sample_id, gsub('^.+_','',run6$sample), sep = "_")
run6$Sample_id = NULL
run6$recoded_Sample_id = NULL


run7$recoded_Sample_id
run7[grepl('Control',run7$Sample_id),]$recoded_Sample_id = c('Dd2_pos1', 'Dd2_pos2')
run7$sample = paste(run7$recoded_Sample_id, gsub('^.+_','',run7$sample), sep = "_")
run7$Sample_id = NULL
run7$recoded_Sample_id = NULL



run8$recoded_Sample_id
run8[grepl('ontrol',run8$Sample_id),]$recoded_Sample_id = c('EXTR-neg1',
                                                             'EXTR-neg2',
                                                             'EXTR-neg3',
                                                             'PCR-neg1',
                                                             'PCR-neg2',
                                                             'Dd2_pos1',
                                                             'Dd2_pos2',
                                                             'Dd2_pos3')
run8 = run8[run8$recoded_Sample_id !="OMIT",]
run8$sample = paste(run8$recoded_Sample_id, gsub('^.+_','',run8$sample), sep = "_")
run8$Sample_id = NULL
run8$recoded_Sample_id = NULL



sum(!(genotyped_samples %in% full_metadata$Sample_id))


genotyped_samples2 = c(run1$sample,
                      run2$sample,
                      run3$sample,
                      run4$sample,
                      run5$sample,
                      run6$sample,
                      run7$sample,
                      run8$sample)

genotyped_samples2 = genotyped_samples2[grepl('ID', ignore.case = T, genotyped_samples2)]
genotyped_samples2 = unique(gsub('_S\\d+$','',genotyped_samples2))


write.table(run1, 'data/sequencing_data/run1/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run2, 'data/sequencing_data/run2/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run3, 'data/sequencing_data/run3/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run4, 'data/sequencing_data/run4/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run5, 'data/sequencing_data/run5/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run6, 'data/sequencing_data/run6/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run7, 'data/sequencing_data/run7/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

write.table(run8, 'data/sequencing_data/run8/dada2/run_dada2/CIGARVariants_Bfilter.out.tsv',
            row.names = F, quote = F)

names(metadata)
metadata = metadata[,c('Sample_id',
                       'Date_of_Collection',
                       'Month_of_Collection',
                       'Quarter_of_Collection',
                       'Year_of_Collection',
                       'Subnational_level0', 
                       'Subnational_level1',
                       'Subnational_level2')]


unique(metadata$Subnational_level2)

metadata[metadata$Subnational_level0 == "Pacific Coast",]$Subnational_level0 = 'Region 1'
metadata[metadata$Subnational_level0 == "East",]$Subnational_level0 = 'Region 2'

metadata[metadata$Subnational_level1 == "Valle del Cauca",]$Subnational_level1 = 'Department 2'
metadata[metadata$Subnational_level1 == "Chocó",]$Subnational_level1 = 'Department 1'
metadata[metadata$Subnational_level1 == "Cauca",]$Subnational_level1 = 'Department 3'
metadata[metadata$Subnational_level1 == "Nariño",]$Subnational_level1 = 'Department 4'
metadata[metadata$Subnational_level1 == "Guainía",]$Subnational_level1 = 'Department 5'

metadata[metadata$Subnational_level2 == "Buenaventura",]$Subnational_level2 = 'Municipality 2'
metadata[metadata$Subnational_level2 == "Quibdó",]$Subnational_level2 = 'Municipality 1'
metadata[metadata$Subnational_level2 == "Guapi",]$Subnational_level2 = 'Municipality 3'
metadata[metadata$Subnational_level2 == "Tumaco",]$Subnational_level2 = 'Municipality 4'
metadata[metadata$Subnational_level2 == "Puerto Inírida",]$Subnational_level2 = 'Municipality 5'


write.csv(metadata, 'data/metadata.csv', row.names = F, quote = F)







