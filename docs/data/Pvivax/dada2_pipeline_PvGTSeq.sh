Rscript ../../../malaria_amplicon_pipeline/create_meta.R -i Fastq -o PvGTSeq_metafile.csv -fw _L001_R1_001.fastq.gz -rv _L001_R2_001.fastq.gz
mkdir dada2
python ../../../malaria_amplicon_pipeline/malaria-amplicon-pipeline/AmpliconPipeline.py \
  --json PvGTSeq_inputs.json
mv prim_fq preprocess_fq run_dada2 prim_meta.txt preprocess_meta.txt stderr.txt stdout.txt ./dada2 