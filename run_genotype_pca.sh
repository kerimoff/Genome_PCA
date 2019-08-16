nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/ImmVar/genotypes/Michigan_GRCh37_1KGPhase3_150119/GRCh38/ImmVar_GRCh38.filtered.vcf.gz\
 --data_name ImmVar\
 -resume
