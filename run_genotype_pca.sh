nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Alasoo_2018/genotypes/Alasoo_2018_GRCh38.filtered.vcf.gz\
 --data_name Alasoo_2018\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Alasoo_2018/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/BLUEPRINT/genotypes/GRCh38/BLUEPRINT_06092016_GRCh38_filtered.vcf.gz\
 --data_name BLUEPRINT\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_BLUEPRINT/\ 
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/BrainSeq/genotypes/Michigan_GRCh37_Phase3_250819/merged/BrainSeq_GRCh38_filtered.vcf.gz\
 --data_name BrainSeq\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_BrainSeq/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
 --data_name Fairfax_2014\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Fairfax_2014/\
 -resume

 nextflow run nextflow_pca_map.nf\
  -profile tartu_hpc\
  --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/GENCORD2/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/GENCORD_GRCh38.filtered.vcf.gz\
  --data_name GENCORD\
  --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_GENCORD/\
  -resume

nextflow run nextflow_pca_map.nf\
  -profile tartu_hpc\
  --vcf /gpfs/hpc/home/a72094/datasets/open_access/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz\
  --data_name GEUVADIS\
  --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_GEUVADIS/\
  -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/HipSci/genotypes/HipSci_GRCh38.filtered.needed.323.samples.vcf.gz\
 --data_name HipSci\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_HipSci_323samples_5PCs/\
 -resume

nextflow run nextflow_pca_map.nf\
  -profile tartu_hpc\
  --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.filtered.vcf.gz\
  --data_name Kasela_2017\
  --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Kasela_2017/\
  -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Lepik_2017/genotypes/GRCh38/Lepik_2017_GRCh38.filtered.vcf.gz\
 --data_name Lepik_2017\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Lepik_2017/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Nedelec_2016/genotypes/Michigan_GRCh37_130618/GRCh38/Nedelec_2016_GRCh38.filtered.vcf.gz\
 --data_name Nedelec_2016\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Nedelec_2016/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Quach_2016/genotypes/Michigan_GRCh37_Phase3_240918/GRCh38/Quach_2016_GRCh38.filtered.vcf.gz\
 --data_name Quach_2016\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Quach_2016/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/Michigan_GRCh37_Phase3_200819/merged/ROSMAP_GRCh38_filtered.vcf.gz\
 --data_name ROSMAP\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_ROSMAP/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Schwartzentruber_2018/Schwartzentruber_2018_needed_samples.vcf.gz\
 --data_name Schwartzentruber_2018\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Schwartzentruber_2018/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Schmiedel_2018/genotypes/Michigan_1KGPhase3_220619/GRCh38/Schmiedel_2018_GRCh38.filtered.vcf.gz\
 --data_name Schmiedel_2018\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_Schmiedel_2018/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/TwinsUK/genotypes/vcf/GRCh38/TwinsUK_GRCh38.filtered.vcf.gz\
 --data_name TwinsUK\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_TwinsUK/\
 -resume

nextflow run nextflow_pca_map.nf\
 -profile tartu_hpc\
 --vcf /gpfs/hpc/home/a72094/datasets/controlled_access/van_de_Bunt_2015/genotypes/Michigan_GRCh37_1KGPhase3_051018/GRCh38/van_de_Bunt_2015_GRCh38.filtered.vcf.gz\
 --data_name van_de_Bunt_2015\
 --outdir /gpfs/hpc/home/kerimov/Genome_QC_Dev/results_van_de_Bunt_2015/\
 -resume

