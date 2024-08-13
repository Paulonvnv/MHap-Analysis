Rscript /Users/pam3650/Documents/Github/malaria_amplicon_pipeline/malaria-amplicon-pipeline/postProc_dada2.R \
  -s dada2/run_dada2/seqtab.tsv \
  --strain PvP01 \
  -ref /Users/pam3650/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq271_refseqs.fasta \
  -b dada2/run_dada2/ASVBimeras.txt \
  -o dada2/run_dada2/ASVTable.txt \
  --fasta \
  --parallel
  
cut -f1,3- dada2/run_dada2/ASVTable.txt > dada2/run_dada2/ASVTable2.txt