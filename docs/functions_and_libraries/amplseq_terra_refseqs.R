# Reference table ----

amplseq_marker_table = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/markers.csv')

amplseq_marker_table %<>% mutate(old_amplicon_name = amplicon,
                                 amplicon = case_when(
                                   amplicon != 'PvDHFR' ~ paste(chromosome, start, end, sep = '_'),
                                   .default = amplicon))

# Reference sequences----

amplseq_ref_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/pf3d7_ref_updated_v3.fasta')

for(old_amplicon in amplseq_marker_table$old_amplicon_name){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name == old_amplicon, ][['amplicon']]
  
  names(amplseq_ref_seq)[names(amplseq_ref_seq) == old_amplicon] = new_name
  
}


sapply(names(amplseq_ref_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
  write(as.character(amplseq_ref_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
})


# Forward primer sequences----

amplseq_fw_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/gtseq_pr_fwd.min_overlap.fasta')

for(old_amplicon in amplseq_marker_table$old_amplicon_name){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name == old_amplicon, ][['amplicon']]
  
  new_name = gsub(paste0('^', old_amplicon, ';'), paste0(new_name, ';'), names(amplseq_fw_seq)[grepl(old_amplicon, names(amplseq_fw_seq))])
  
  names(amplseq_fw_seq)[grepl(old_amplicon, names(amplseq_fw_seq))] = new_name
  
}

#names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))] = gsub('_NA_NA', 'PvDHFR', names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))])

sapply(names(amplseq_fw_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_fw_seqs.fasta', append = T)
  write(as.character(amplseq_fw_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_fw_seqs.fasta', append = T)
})

# Reverse primer sequences----

amplseq_rv_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/gtseq_pr_fwd.min_overlap.fasta')

for(old_amplicon in amplseq_marker_table$old_amplicon_name){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name == old_amplicon, ][['amplicon']]
  
  new_name = gsub(paste0('^', old_amplicon, ';'), paste0(new_name, ';'), names(amplseq_rv_seq)[grepl(old_amplicon, names(amplseq_rv_seq))])
  
  names(amplseq_rv_seq)[grepl(old_amplicon, names(amplseq_rv_seq))] = new_name
  
}

#names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))] = gsub('_NA_NA', 'PvDHFR', names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))])

sapply(names(amplseq_rv_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_rv_seqs.fasta', append = T)
  write(as.character(amplseq_rv_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_rv_seqs.fasta', append = T)
})



