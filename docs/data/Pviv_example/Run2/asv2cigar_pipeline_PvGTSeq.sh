python /Users/pam3650/Documents/Github/malaria_amplicon_pipeline/malaria-amplicon-pipeline/ASV_to_CIGAR.py \
  -d /Users/pam3650/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_refseqs.fasta \
  --asv_to_cigar dada2/run_dada2/ASV_to_CIGAR.out.txt \
  dada2/run_dada2/ASVSeqs.fasta \
  dada2/run_dada2/ASVTable2.txt \
  dada2/run_dada2/seqtab.tsv \
  dada2/run_dada2/CIGARVariants_Bfilter.out.tsv \
  --exclude_bimeras