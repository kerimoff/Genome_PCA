nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/BrainSeq/genotypes/Michigan_GRCh37_Phase3_250819/merged/BrainSeq_GRCh38_filtered.vcf.gz\
 --data_name BrainSeq\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_BrainSeq/\
 -resume

