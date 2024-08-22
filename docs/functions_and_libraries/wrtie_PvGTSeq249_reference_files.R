

pvgtseq_markers271 = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_markersTable.csv')
pvgtseq_fw271 = Biostrings::readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_fwd.fasta')
pvgtseq_rv271 = Biostrings::readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_rvs.fasta')
pvgtseq_refseqs271 = Biostrings::readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_refseqs.fasta')

pvgtseq_markers249 = pvgtseq_markers271 %>% filter(grepl('249',Set))

write.csv(pvgtseq_markers249, 
          '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_markersTable.csv',
          quote = F,
          row.names = F)


pvgtseq_fw249 = pvgtseq_fw271[grepl(paste0('(', paste(pvgtseq_markers249$amplicon, collapse = ';|'), ';)'),
                                    names(pvgtseq_fw271))]

sapply(names(pvgtseq_fw249), function(amplicon){
  write(paste0('>', amplicon),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_fwd.fasta', 
        sep = '\t', 
        append = TRUE)
  
  write(paste0('^',as.character(pvgtseq_fw249[[amplicon]])),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_fwd.fasta', 
        sep = '\t', 
        append = TRUE)
})

pvgtseq_rv249 = pvgtseq_rv271[grepl(paste0('(', paste(pvgtseq_markers249$amplicon, collapse = ';|'), ';)'),
                                    names(pvgtseq_rv271))]

sapply(names(pvgtseq_rv249), function(amplicon){
  write(paste0('>', amplicon),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_rvs.fasta', 
        sep = '\t', 
        append = TRUE)
  
  write(paste0('^',as.character(pvgtseq_rv249[[amplicon]])),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_rvs.fasta', 
        sep = '\t', 
        append = TRUE)
})


pvgtseq_refseqs249 = pvgtseq_refseqs271[grepl(paste0('(^', paste(pvgtseq_markers249$amplicon, collapse = '$|^'), '$)'),
                                    names(pvgtseq_refseqs271))]

sapply(names(pvgtseq_refseqs249), function(amplicon){
  write(paste0('>', amplicon),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_refseqs.fasta', 
        sep = '\t', 
        append = TRUE)
  
  write(as.character(pvgtseq_refseqs249[[amplicon]]),
        file = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_refseqs.fasta', 
        sep = '\t', 
        append = TRUE)
})

