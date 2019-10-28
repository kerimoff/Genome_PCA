#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
options(bitmapType='cairo')

library(ggplot2)
library(dplyr)
library(DMwR)
library("data.table")

generate_distance_matrix <- function(reference_pca_df, pca_df, n_pcs = 3, method = "euclidean"){
  # required fields assertion ====
  assertthat::assert_that(hasName(pca_df, "genotype_id"), msg = "Column genotype_id is missing in pca_df")
  assertthat::assert_that(hasName(pca_df, "superpopulation_code"), msg = "Column superpopulation_code is missing in pca_df")
  pcs <- c()
  for (i in 1:n_pcs) {
    pc_num <- paste0("PC",i)
    assertthat::assert_that(hasName(pca_df, pc_num), msg = paste0("Column ", pc_num, " is missing in pca_df"))
    pcs <- c(pcs, pc_num)
  }
  
  assertthat::assert_that(hasName(reference_pca_df, "genotype_id"), msg = "Column genotype_id is missing in reference_pca_df")
  assertthat::assert_that(hasName(reference_pca_df, "superpopulation_code"), msg = "Column superpopulation_code is missing in reference_pca_df")
  for (i in 1:n_pcs) {
    pc_num <- paste0("PC",i)
    assertthat::assert_that(hasName(reference_pca_df, pc_num), msg = paste0("Column ", pc_num, " is missing in reference_pca_df"))
  }
  # ====
  
  sample_size <- pca_df %>% nrow()
  reference_pca_df$genotype_id <- paste0("ref_", reference_pca_df$genotype_id)
  comb_pcs = rbind(pca_df, reference_pca_df)
  rownames(comb_pcs) <- comb_pcs$genotype_id
  distance_matrix <- dist(comb_pcs[pcs], method = method) %>% as.matrix(labels = TRUE)
  distance_matrix <- distance_matrix[-c(1:sample_size), 1:sample_size] 
  tibble_ind <- distance_matrix %>% as_tibble(rownames="genotype_id")
  tibble_ind <- reference_pca_df[,c("genotype_id","superpopulation_code")] %>% left_join(tibble_ind)
  
  mean_ind_diff_summary <- tibble_ind[,-1] %>% group_by(superpopulation_code) %>% summarize_all("mean")
  return(mean_ind_diff_summary)
}

assign_populations <- function(distance_matrix, abs_threshold = 0.02, rel_threshold = 1.7){
  assigned <- distance_matrix[,-1] %>% 
    rbind(if_else(apply(distance_matrix[,-1], 2, min)<abs_threshold, 
                  as.character(distance_matrix$superpopulation_code[apply(distance_matrix[,-1], 2, which.min)]), 
                  "Admixed")) %>% as.data.frame()
  rownames(assigned) <- c(as.character(distance_matrix$superpopulation_code), "pop_assign_abs_thresh")
  assigned <- assigned %>% t() %>% as.data.frame() 
  assigned <- cbind(genotype_id=rownames(assigned), assigned)

  # find if there are admixed samples
  mins <- apply(distance_matrix[,-1], 2, sort) %>% t() %>% as.data.frame() 
  mins <- mins %>% mutate(isAdmixed = mins[,1]*rel_threshold > mins[,2])
  # set assigned_population as minimum distance population first
  assigned_rel <- distance_matrix[,-1] %>% rbind(as.character(distance_matrix$superpopulation_code[apply(distance_matrix[,-1], 2, which.min)])) %>% as.data.frame()
  rownames(assigned_rel) <- c(as.character(distance_matrix$superpopulation_code), "pop_assign_rel_thresh")
  assigned_rel <- assigned_rel %>% t() %>% as.data.frame()
  
  # if there are admixed samples according to proportional threshold replace assigned_rel population value
  if (sum(mins$isAdmixed)>0) {
    levels(assigned_rel$pop_assign_rel_thresh) <- c(levels(assigned_rel$pop_assign_rel_thresh), "Admixed")
    assigned_rel$pop_assign_rel_thresh[mins$isAdmixed] <- "Admixed"
  }
  assigned_rel <- cbind(genotype_id=rownames(assigned_rel), assigned_rel)

  return(left_join(assigned, assigned_rel[,c("genotype_id", "pop_assign_rel_thresh")]))
}

if (FALSE) {
  source_populations_file <- "temp/fairfax/igsr_samples.tsv"
  pca_table <- "temp/fairfax/main_overlapped_pca.vect"
  projections <- "temp/fairfax/new_dataset_scores.profile.adj"
  data_name <- "fairfax"
  distance_method <- "euclidean"
  n_pcs <- 5
  admixed_abs_threshold <- 0.02
  admixed_rel_threshold <- 1.7
}

pca_table = if (is.na(args[1])) 'main_overlapped_pca.vect' else args[1]
projections = if (is.na(args[2])) 'new_dataset_scores.profile.adj' else args[2]
source_populations_file = if (is.na(args[3])) "samples_data.tsv" else args[3]
data_name = if (is.na(args[4])) 'New_dataset' else args[4] 
n_pcs = if (is.na(args[5])) 3 else args[5] 
distance_method = if (is.na(args[6])) "euclidean" else args[6] 
admixed_abs_threshold = if (is.na(args[7])) 0.02 else args[7] 
admixed_rel_threshold = if (is.na(args[8])) 1.7 else args[8] 

if (!dir.exists("plots")) { dir.create("plots") }

data_colnames <- c("genotype_id")
for (i in 1:n_pcs) {
  data_colnames <- c(data_colnames, paste0("PC",i))
}
pca <- read.table(pca_table, header = FALSE, )
pca <- pca %>% select(-c(2)) %>% setNames(data_colnames)

source_populations = read.table(source_populations_file, header = TRUE, sep='\t')
main_pca <- merge(pca, select(source_populations, genotype_id, superpopulation_code), by  = "genotype_id") 

ggplot(main_pca, aes(x=PC1, y=PC2, color=as.factor(superpopulation_code))) +
  geom_point() + guides(color = guide_legend(title='Superpopulation')) + 
  labs(x = "PC1", y = "PC2") + ggtitle(data_name)+ coord_fixed()

ggsave('main_pca.pdf', path = "plots", device = "pdf")

projections_pcs <- read.table(projections, header = TRUE)
if (sum(projections_pcs$ID1==projections_pcs$ID2)==0) {
  projections_pcs$ID1 <- paste0(projections_pcs$ID1, "_", projections_pcs$ID2)
}
projections_pcs <- projections_pcs %>% select(-c(2,3,4)) %>% 
  setNames(data_colnames) %>% 
  mutate(superpopulation_code = data_name)

ggplot(projections_pcs, aes(x=PC1, y=PC2)) + geom_point()+ coord_fixed()
ggsave('projections_only.pdf', path = "plots", device = "pdf")

comb_pcs = rbind(main_pca, projections_pcs)

ggplot(comb_pcs, aes(x=PC1, y=PC2, color=as.factor(superpopulation_code))) + geom_point() + 
  guides(color = guide_legend(title='Superpopulation')) + ggtitle(data_name) + coord_fixed() 
ggsave('projections_on_main.pdf', path = "plots", device = "pdf")

################### Populations assignment ################### 
means_main_pca <- main_pca %>% group_by(superpopulation_code) %>% summarize(mean_P1 = mean(PC1), mean_P2 = mean(PC2))

distance <- function(p1, p2, mean1, mean2){
  return(sqrt((mean1-p1)**2 + (mean2-p2)**2))
}

populations = unique(main_pca$superpopulation_code)
samples_with_dists = projections_pcs

for (population_code in populations) {
  col_name = sprintf('dist_%s', population_code)
  samples_with_dists = samples_with_dists %>% 
    mutate( !!col_name := distance(PC1, PC2, means_main_pca[means_main_pca$superpopulation_code == population_code, ]$mean_P1,
                                   means_main_pca[means_main_pca$superpopulation_code == population_code, ]$mean_P2))
}

# samples_with_dists %>% mutate(sum = select(., cols) %>% apply(1, sum))
samples_with_dists = samples_with_dists %>% mutate(sum_dist = select(., starts_with('dist_')) %>% apply(1, sum))

for (population_code in populations) {
  col_name = sprintf('dist_%s', population_code)
  col_name_norm = sprintf('dist_%s_norm', population_code)
  samples_with_dists = samples_with_dists %>%
    mutate( !!col_name_norm := !!as.name(col_name) / sum_dist)
}
############### k-nearest neighbours #############

train <- main_pca %>% select(PC1, PC2, superpopulation_code)
test <- samples_with_dists %>% select(PC1, PC2)

knn_pops <- kNN(superpopulation_code ~ PC1+PC2,train,test,norm=FALSE,k=5)
new_populations <- samples_with_dists %>% mutate(knn_pop = knn_pops)

main_pca_to_map <- main_pca %>% mutate(knn_pop = superpopulation_code)

ggplot(new_populations, aes(x=PC1, y=PC2, color=as.factor(knn_pop))) + 
  geom_point() + 
  geom_point(data=main_pca_to_map, alpha = 0.3, shape=4) + 
  guides(color = guide_legend(title='KNN population')) + ggtitle(data_name) + coord_fixed() 
ggsave('knn.pdf', path = "plots", device = "pdf")

threshold = 0.1

new_populations_threshold = new_populations %>% 
  mutate(min_val = apply(select(new_populations, ends_with('norm')), 1, min)) %>% 
  mutate(knn_pop = if_else(min_val>threshold, 'АdMixed', as.character(knn_pop)))  %>% 
  select(-min_val)

ggplot(new_populations_threshold, aes(x=PC1, y=PC2, color=as.factor(knn_pop))) + 
  geom_point() + 
  geom_point(data=main_pca_to_map, alpha = 0.3, shape=4) + 
  guides(color = guide_legend(title='KNN population')) + ggtitle(data_name) + coord_fixed() 

ggsave('knn_threshold.pdf', path = "plots", device = "pdf")

###################################################
finalData <- new_populations %>% mutate(knn_pop_threshold = new_populations_threshold$knn_pop)
write.table(finalData, file='populations.tsv', quote = FALSE, row.names = FALSE, sep = '\t')

############### distance matrix generation and populain assignment with threshold  #############

for (pc_index in 2:5) {
  distance_matrix <- generate_distance_matrix(reference_pca_df = main_pca, pca_df = projections_pcs, method = distance_method, n_pcs = pc_index)  
  summ_diff_matrix_melted <- melt(distance_matrix)
  dist_plot <- ggplot(summ_diff_matrix_melted, aes(value, fill = superpopulation_code)) + geom_density(alpha = 0.2)
  ggsave(filename = paste0("distance_histogram_", pc_index, 'PCs.pdf'), plot = dist_plot, path = "plots", device = "pdf")
}

distance_matrix <- generate_distance_matrix(reference_pca_df = main_pca, pca_df = projections_pcs, method = distance_method, n_pcs = 3)

assigned_populations <- assign_populations(distance_matrix, abs_threshold = admixed_abs_threshold, rel_threshold = admixed_rel_threshold)
filename <- paste0("pop_assigned_abs_", admixed_abs_threshold, "_rel_", admixed_rel_threshold, ".tsv")
write.table(assigned_populations, file=filename, quote = FALSE, row.names = FALSE, sep = '\t', col.names = TRUE)

main_pca %>% mutate()
assigned_with_pcs <- projections_pcs %>% select(-c("superpopulation_code")) %>% left_join(assigned_populations[c("genotype_id","pop_assign_abs_thresh")])
colnames(assigned_with_pcs)[which(colnames(assigned_with_pcs)=="pop_assign_abs_thresh")] <- "superpopulation_code"
assigned_populations_abs_tresh <- ggplot(assigned_with_pcs, aes(x=PC1, y=PC2, color=superpopulation_code)) + 
  geom_point() + 
  geom_point(data=main_pca, alpha = 0.1, shape=4) + 
  guides(color = guide_legend(title='Assigned population\n(abs threshold)')) + ggtitle(data_name) + coord_fixed() 
ggsave(filename = paste0(data_name, "_pop_assign_abs_treshold_", admixed_abs_threshold, ".pdf"), plot = assigned_populations_abs_tresh, path = "plots", device = "pdf")

assigned_with_pcs <- projections_pcs %>% select(-c("superpopulation_code")) %>% left_join(assigned_populations[c("genotype_id","pop_assign_rel_thresh")])
colnames(assigned_with_pcs)[which(colnames(assigned_with_pcs)=="pop_assign_rel_thresh")] <- "superpopulation_code"
assigned_populations_rel_tresh <- ggplot(assigned_with_pcs, aes(x=PC1, y=PC2, color=superpopulation_code)) + 
  geom_point() + 
  geom_point(data=main_pca, alpha = 0.1, shape=4) + 
  guides(color = guide_legend(title='Assigned population\n(rel threshold)')) + ggtitle(data_name) + coord_fixed() 
ggsave(filename = paste0(data_name,"_pop_assign_rel_treshold_", admixed_rel_threshold, ".pdf"), plot = assigned_populations_rel_tresh, path = "plots", device = "pdf")


