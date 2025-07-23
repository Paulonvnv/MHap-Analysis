

MAD4HAtTer_homopolymers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/MAD4HAtTer_homopolymers.csv')
ref_fasta = "~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/MAD4HatTer_refseqs_v4.fasta"
homopolymer_length = 5

mhaps = MAD4HAtTer_homopolymers



mhaps$homopolymer_regions = NA
mhaps$homopolymer_anchors = NA

ref_sequences = readDNAStringSet(ref_fasta)

if(sum(grepl('(/|-|:)', names(ref_sequences))) > 0){
  
  names(ref_sequences) = gsub('(/|-|:)', '_', names(ref_sequences))
  print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the names of the fasta sequences.')
  
}

homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'

homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)


for(mhap in mhaps$amplicon){
  
  # Identification of Homopolymers
  
  homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
  
  if(length(homopolymers) > 0){
    
    homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
    
    mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                           paste(homopolymers_location[[1]][,'start'],
                                                                                 homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                     collapse = ',')
    
    homopolymer_starts = homopolymers_location[[1]][,'start']
    
    anchor_position = NULL
    for(homopolymer in 1:length(homopolymers_location[[1]][,'start'])){
      
      homopolymer_start = homopolymers_location[[1]][,'start'][homopolymer]
      
      if(homopolymer_start > 1){
        
        anchor_position = c(anchor_position, 
                            as.character(Biostrings::substr(ref_sequences[[mhap]], 
                                                            start = homopolymer_start - 1, 
                                                            stop = homopolymer_start - 1)))
      }else{
        anchor_position = c(anchor_position, '')
      }
      
    }
    
    if(!is.null(anchor_position)){
      homopolymer_anchors = paste0(homopolymers_location[[1]][,'start'] - 1, anchor_position)
      mhaps[mhaps$amplicon == mhap, ][['homopolymer_anchors']] = paste(homopolymer_anchors, collapse = ',')  
    }
    
  }
  
}

MAD4HAtTer_homopolymers$homopolymer_regions = mhaps$homopolymer_regions

MAD4HAtTer_homopolymers$homopolymer_anchors = mhaps$homopolymer_anchors


write.csv(MAD4HAtTer_homopolymers, '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/MAD4HAtTer_homopolymers.csv', quote = F, row.names = F)


