#!/usr/bin/env nextflow

// executables 
plink = "/gpfs/hpchome/peikova/software/plink/plink"
ldak5 = "/gpfs/hpchome/peikova/software/ldak5.linux"
// path to script for plotting and mapping
plot_pca = "plot_pca.R"

params.data_name = 'ImmVar'
params.vcf = "/gpfs/hpchome/a72094/hpc/projects/population_PCA/genotypes/ImmVar_GRCh38.filtered.vcf.gz"

// file with populations of main dataset
params.populations_file = "/gpfs/hpchome/peikova/bioinf2/project/igsr_samples.tsv"

// path to main dataset
params.main_vcf = "/gpfs/hpchome/peikova/bioinf2/project/source_data/GRCh38_renamed_ids_no_multiallelic.vcf.gz"

params.num_pc = 3
params.threads = 16
params.outdir = "nf_results/$params.data_name"

// file with samples (americans) to remove for pca mapping
params.ids_to_remove_file = "/gpfs/hpchome/peikova/bioinf2/project/amrs.txt"
// flag if remove samples before mapping
params.exclude_population = true

vcf_file = file(params.vcf)
main_vcf = file(params.main_vcf)
populations_file = file(params.populations_file)

// TODO:    
// manage arguments and used programs
// add qtltools option
// add option for only pca without mapping

process main_vcf_to_binary{
    input:
    file "vcf_main.vcf.gz" from main_vcf

    output:
    set file('main.bed'), file('main.bim'), file('main.fam') into main_to_delete_dublicates

    script:
    """
    # do ld pruning and save resulted snps in file  
    $plink --vcf vcf_main.vcf.gz --vcf-half-call h --indep-pairwise 50000 200 0.05 --out main_pruned_varaints_list --threads $params.threads

    # make bfiles for pruned 1000 genome proj 
    $plink --vcf vcf_main.vcf.gz --vcf-half-call h --extract main_pruned_varaints_list.prune.in --make-bed --out main
    """
}

if(params.exclude_population){
    process remove_family{
    input:
    set file('main.bed'), file('main.bim'), file('main.fam') from main_to_delete_dublicates
    file 'ids_to_remove.txt' from file(params.ids_to_remove_file)

    output:
    set file('main_no_dubl.bed'), file('main_no_dubl.bim'), file('main_no_dubl.fam') into  main_to_extract_snps, main_binary_source

    script:
    """
    $plink --bfile main --remove-fam ids_to_remove.txt --make-bed --out main
    # finds dublicate vars
    $plink --bfile main --list-duplicate-vars --out dubl

    # delete dublicate vars
    $plink --bfile main --exclude dubl.dupvar --snps-only --make-bed --out main_no_dubl
    """
    }
}else{
    process remove_dubl{
    input:
    set file('main.bed'), file('main.bim'), file('main.fam') from main_to_delete_dublicates
  
    output:
    set file('main_no_dubl.bed'), file('main_no_dubl.bim'), file('main_no_dubl.fam') into main_to_extract_snps, main_binary_source

    script:
    """
    # finds dublicate vars
    $plink --bfile main --list-duplicate-vars --out dubl
    # delete dublicate vars
    $plink --bfile main --exclude dubl.dupvar --snps-only --make-bed --out main_no_dubl
    """
    }

}

// convert vcf file to plink binary file (.bed)
process convertVCFtoBED{
    input:
    file 'source.vcf.gz' from vcf_file

    output:
    set file ('binary_source.bed'), file('binary_source.bim'), file('binary_source.fam') into bed_files

    script:
    """
    $plink --vcf source.vcf.gz --out binary_source --threads $params.threads
    
    $plink --bfile binary_source --list-duplicate-vars --out list_dubl
    $plink --bfile binary_source --exclude list_dubl.dupvar --snps-only --make-bed --out binary_source
    """
}


process get_SNPs_list_from_main_dataset{
    input: 
    set file('main_source.bed'), file('main_source.bim'), file('main_source.fam') from main_to_extract_snps

    output:
    file 'main_snps_list.snplist' into main_snps_file

    script:
    """
    $plink --bfile main_source --write-snplist --out main_snps_list --snps-only
    """
}

process extract_SNPs_and_make_bed{
    input:
    set file ('binary_source.bed'), file('binary_source.bim'), file('binary_source.fam') from bed_files
    file 'main.snplist' from main_snps_file

    output:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') into new_dataset_overlapped1, new_dataset_overlapped2
    file 'overlapped_snps.snplist' into new_dataset_overlapped_snplist

    // extract snps present in pruned data from new dataset
    // --make-bed makes sure bfiles created!
    """
    $plink --bfile binary_source --extract main.snplist --make-bed --out new_dataset_overlapped
    $plink --bfile new_dataset_overlapped --write-snplist --out overlapped_snps
    """
}

process extract_overlapped_SNPs_from_main_dataset {
    input:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') from new_dataset_overlapped1
    file 'overlapped_snps.snplist' from new_dataset_overlapped_snplist
    set file('main_source.bed'), file('main_source.bim'), file('main_source.fam') from main_binary_source

    output:
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') into main_to_kins
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') into main_bed_to_pca

    """
    $plink --bfile main_source --extract overlapped_snps.snplist --make-bed --out main_overlapped
    """
}

process calc_kins_matrices{
    input:
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') from main_to_kins

    output:
    set file("main_overlapped_kins.grm.bin"), file("main_overlapped_kins.grm.id"), file("main_overlapped_kins.grm.adjust"), file("main_overlapped_kins.grm.details") into main_to_pca
    
    """
    $ldak5  --calc-kins-direct main_overlapped_kins --bfile main_overlapped --ignore-weights YES --power -0.25
    """
}

process calc_pca_and_loads{
    publishDir "$params.outdir", mode: 'copy'

    input:
    set file("main_overlapped_kins.grm.bin"), file("main_overlapped_kins.grm.id"), file("main_overlapped_kins.grm.adjust"), file("main_overlapped_kins.grm.details") from main_to_pca
    set file('main_overlapped.bed'), file('main_overlapped.bim'), file('main_overlapped.fam') from main_bed_to_pca

    output:
    file 'main_overlapped_loads.load' into loads_for_mapping
    file 'main_overlapped_pca.vect' into plot_data_vect
     
    """
    $ldak5 --pca main_overlapped_pca --grm main_overlapped_kins --axes $params.num_pc
    $ldak5 --calc-pca-loads main_overlapped_loads --grm main_overlapped_kins --pcastem main_overlapped_pca --bfile main_overlapped
    """
}

process map_new_dataset{
    publishDir "$params.outdir", mode: 'copy'
    input:
    set file('new_dataset_overlapped.bed'), file('new_dataset_overlapped.bim'), file('new_dataset_overlapped.fam') from new_dataset_overlapped2
    file 'main_overlapped_loads.load' from loads_for_mapping

    output:
    file 'new_dataset_scores.profile.adj' into plot_data

    """
    $ldak5  --calc-scores new_dataset_scores --bfile new_dataset_overlapped --scorefile main_overlapped_loads.load --power 0
    """

}

process plot_pca{
    publishDir "$params.outdir", mode: 'copy'
    input:
    file 'new_dataset_scores.profile.adj' from plot_data
    file 'main_overlapped_pca.vect' from plot_data_vect
    file 'samples_data.tsv' from populations_file

    output:
    set file('main_pca.png'), file('projections_only.png'), file('projections_on_main.png'), file('populations.tsv'), file('knn_threshold.png'), file('knn.png')

    """
    Rscript $plot_pca main_overlapped_pca.vect new_dataset_scores.profile.adj samples_data.tsv $params.data_name
    """

}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}






