
library(magrittr)
library(tidyr)
library(dplyr)

pfmarkers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/markers.csv')


AMPLSeq_coordinates = pfmarkers %>% 
  filter(amplicon != 'PvDHFR') %>%
  select(amplicon, chromosome, start, end) %>%
  mutate(start = as.integer(start),
         end = as.integer(end))

AMPLSeq_positions = NULL

for(amplicon in AMPLSeq_coordinates$amplicon){
  
  AMPLSeq_positions = rbind(
    AMPLSeq_positions,
    data.frame(CHROM = AMPLSeq_coordinates[AMPLSeq_coordinates$amplicon == amplicon, ][['chromosome']],
               POS = AMPLSeq_coordinates[AMPLSeq_coordinates$amplicon == amplicon, ][['start']]:AMPLSeq_coordinates[AMPLSeq_coordinates$amplicon == amplicon, ][['end']]
               ))
  }




AMPLSeq_positions %<>% rename('CHROM' = 'chromosome', 'POS' = 'pos')

write.table(AMPLSeq_positions, 
            '~/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/AMPLSeq_positions.bed', 
            sep = '\t', 
            quote = F,
            row.names = F)
