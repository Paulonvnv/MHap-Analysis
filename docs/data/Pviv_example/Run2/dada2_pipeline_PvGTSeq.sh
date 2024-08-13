Rscript /Users/pam3650/Documents/Github/malaria_amplicon_pipeline/copy_MiSeq_Fastqs.R -i Pv_Run1-413086072/FASTQ_Generation_2024-03-27_10_13_35Z-725582948 -o Fastq
Rscript /Users/pam3650/Documents/Github/malaria_amplicon_pipeline/create_meta.R -i Fastq -o PvGTSeq_metafile.csv -fw _L001_R1_001.fastq.gz -rv _L001_R2_001.fastq.gz
mkdir dada2
python /Users/pam3650/Documents/Github/malaria_amplicon_pipeline/malaria-amplicon-pipeline/AmpliconPipeline.py \
  --json PvGTSeq_inputs.json
mv prim_fq preprocess_fq run_dada2 prim_meta.txt preprocess_meta.txt stderr.txt stdout.txt ./dada2 