amplseq_ref_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/pf3d7_ref_updated_v3.fasta')

amplseq_marker_table = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/markers.csv')

amplseq_marker_table %<>% mutate(old_amplicon_name = amplicon,
                                 amplicon = paste(chromosome, start, end, sep = '_'))

old_amplicon = amplseq_marker_table$old_amplicon_name[1]
for(old_amplicon in amplseq_marker_table$old_amplicon_name){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name == old_amplicon, ][['amplicon']]
  
  names(amplseq_ref_seq)[names(amplseq_ref_seq) == old_amplicon] = new_name
  
}

names(amplseq_ref_seq)[names(amplseq_ref_seq) == '_NA_NA'] = 'PvDHFR'

sapply(names(amplseq_ref_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
  write(as.character(amplseq_ref_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
})
