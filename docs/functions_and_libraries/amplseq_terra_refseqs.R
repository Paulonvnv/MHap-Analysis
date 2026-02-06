# Reference table ----

amplseq_marker_table = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/markers.csv')

amplseq_marker_table %<>% mutate(old_amplicon_name = amplicon,
                                 amplicon = case_when(
                                   amplicon != 'PvDHFR' ~ paste(chromosome, start, end, sep = '_'),
                                   .default = amplicon))

# Upload old reference files----
amplseq_ref_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/pf3d7_ref_updated_v3.fasta')

amplseq_fw_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/gtseq_pr_fwd.min_overlap.fasta')

amplseq_rv_seq = readDNAStringSet('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/gtseq_pr_rvs.min_overlap.fasta')


not_in_refseqs = names(amplseq_ref_seq)[!(names(amplseq_ref_seq) %in% amplseq_marker_table$old_amplicon_name)]
not_in_marker_table_from_refseqs = amplseq_marker_table$old_amplicon_name[!(amplseq_marker_table$old_amplicon_name %in% names(amplseq_ref_seq))]

not_in_marker_table = gsub(';.+$', '', names(amplseq_fw_seq))[!(gsub(';.+$', '', names(amplseq_fw_seq)) %in% amplseq_marker_table$old_amplicon_name)]
not_in_fwd_seq = amplseq_marker_table$old_amplicon_name[!(amplseq_marker_table$old_amplicon_name %in% gsub(';.+$', '', names(amplseq_fw_seq)))]

not_in_marker_table_from_rev = gsub(';.+$', '', names(amplseq_rv_seq))[!(gsub(';.+$', '', names(amplseq_rv_seq)) %in% amplseq_marker_table$old_amplicon_name_primers)]
not_in_rev_seq = amplseq_marker_table$old_amplicon_name_primers[!(amplseq_marker_table$old_amplicon_name_primers %in% gsub(';.+$', '', names(amplseq_rv_seq)))]

amplseq_marker_table %<>% 
  mutate(old_amplicon_name_primers_fwd = 
           case_when(
             old_amplicon_name == 'PfSERA3' ~ 'SERA3',
             old_amplicon_name == 'PfSERA2_1' ~ 'PfSERA2',
             old_amplicon_name == 'PfCSP_1' ~ 'PfCSP',
             old_amplicon_name == 'PfDHFR_1' ~ 'PfDHFR_16_59',
             old_amplicon_name == 'PfDHFR_2' ~ 'PfDHFR_108_164',
             old_amplicon_name == 'PfMDR1_1' ~ 'MDR1_86',
             old_amplicon_name == 'PfMDR1_2' ~ 'MDR1_184',
             old_amplicon_name == 'PfMDR1_3' ~ 'MDR1_1034_1042',
             old_amplicon_name == 'PfMDR1_4' ~ 'MDR1_1246',
             old_amplicon_name == 'PfDHPS_2' ~ 'PfDHPS_540_613',
             old_amplicon_name == 'PfDHPS_1' ~ 'DHPS_436_437',
             old_amplicon_name == 'PfAMA_1' ~ 'PfAMA1',
             old_amplicon_name == 'PfTRAP_1' ~ 'PfTRAP',
             old_amplicon_name == 'PfKelch13C580Y' ~ 'PfKelch13_520_580',
             .default = old_amplicon_name
             ),
         old_amplicon_name_primers_rev = case_when(
           old_amplicon_name_primers_fwd == 'PfSERA2' ~ 'PfTSERA2',
           old_amplicon_name_primers_fwd == 'MDR1_86' ~ 'MDR1_86_',
           .default = old_amplicon_name_primers_fwd
         )
         )

# Reference sequences----

for(old_amplicon in amplseq_marker_table$old_amplicon_name){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name == old_amplicon, ][['amplicon']]
  
  names(amplseq_ref_seq)[names(amplseq_ref_seq) == old_amplicon] = new_name
  
}


sapply(names(amplseq_ref_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
  write(as.character(amplseq_ref_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_ref_seqs.fasta', append = T)
})


# Forward primer sequences----

for(old_amplicon in amplseq_marker_table$old_amplicon_name_primers_fwd){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name_primers_fwd == old_amplicon, ][['amplicon']]
  
  new_name = gsub(paste0('^', old_amplicon, ';'), paste0(new_name, ';'), names(amplseq_fw_seq)[grepl(old_amplicon, names(amplseq_fw_seq))])
  
  names(amplseq_fw_seq)[grepl(old_amplicon, names(amplseq_fw_seq))] = new_name
  
}

#names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))] = gsub('_NA_NA', 'PvDHFR', names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))])

sapply(names(amplseq_fw_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_fw_seqs.fasta', append = T)
  write(as.character(amplseq_fw_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_fw_seqs.fasta', append = T)
})

# Reverse primer sequences----



for(old_amplicon in amplseq_marker_table$old_amplicon_name_primers_rev){
  
  new_name = amplseq_marker_table[amplseq_marker_table$old_amplicon_name_primers_rev == old_amplicon, ][['amplicon']]
  
  new_name = gsub(paste0('^', old_amplicon, ';'), paste0(new_name, ';'), names(amplseq_rv_seq)[grepl(old_amplicon, names(amplseq_rv_seq))])
  
  names(amplseq_rv_seq)[grepl(old_amplicon, names(amplseq_rv_seq))] = new_name
  
}

#names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))] = gsub('_NA_NA', 'PvDHFR', names(amplseq_fw_seq)[grepl('_NA_NA', names(amplseq_fw_seq))])

sapply(names(amplseq_rv_seq), function(amplicon){
  write(paste0('>', amplicon), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_rv_seqs.fasta', append = T)
  write(as.character(amplseq_rv_seq[[amplicon]]), file = '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_rv_seqs.fasta', append = T)
})



write.csv(amplseq_marker_table, '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/amplseq_marker_table.csv', quote = F, row.names = F)

