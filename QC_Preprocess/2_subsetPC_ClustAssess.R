#First load gtf and get pc IDs
library(dplyr)
library(Seurat)
library(ClustAssess, lib.loc = "/usr/local/lib/R/site-library")
# I used this to get genes
# genes.pc.gtf obrained from gtf by grep for gene_type "protein_coding"
#cat genes.pc.gtf |grep -v '^#'|grep -P '^\S+\t\S+\tgene\S*\t' | cut -d ';' -f4 | cut -d' ' -f3 | sed 's/^"//;s/"$//'
makeCA <- function(so,ncores,outdir_obj,outdir_app,assay_name='SCT',project_name){
  print(project_name)
  so <- SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)
  saveRDS(so,paste0(outdir_obj,project_name,'-so.rds'))
  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )
  
  features <- dimnames(so@assays[[assay_name]])[[1]]
  var_features <- so@assays[[assay_name]]@var.features
  n_abundant <- 3000
  most_abundant_genes <- rownames(so@assays[[assay_name]])[order(Matrix::rowSums(so@assays[[assay_name]]),
                                                                 decreasing = TRUE
  )][1:n_abundant]
  
  gene_list = list(
    "Highly_Variable" = so@assays[[assay_name]]@var.features,
    "Most_Abundant" = most_abundant_genes
  )
  
  steps_list = list(
    "Highly_Variable" = c(500,750,1000,1250,1500,1750,2000),
    "Most_Abundant" =c(500,750,1000,1250,1500,1750,2000)
  )
  
  expr_matrix = so@assays$SCT@scale.data
  gc()
  
  doParallel::registerDoParallel(cl = my_cluster)
  
  RhpcBLASctl::blas_set_num_threads(ncores)
  clustassess_autom <- automatic_stability_assessment(
    expression_matrix = expr_matrix,
    n_repetitions = 100,
    temp_file = 'mytemp_diff2.rds',
    n_neigh_sequence = seq(from = 5, to = 50, by = 5),
    resolution_sequence = seq(from = 0.1, to = 2, by = 0.1),
    features_sets = gene_list,
    steps = steps_list,
    n_top_configs = 7,
    npcs = 30,
    save_temp = TRUE,
    verbose = TRUE,
    umap_arguments = list(
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine",
      init = "spectral"
    )
  )
  
  
  
  saveRDS(clustassess_autom, paste0(outdir_obj,project_name,'-ca.rds')) 
  
  parallel::stopCluster(cl = my_cluster)
  
  ClustAssess::write_shiny_app(
    object = so,
    assay_name = "SCT",
    clustassess_object = clustassess_autom,
    project_folder = paste0(outdir_app,project_name,'/'),
    shiny_app_title = paste0('Lottie - ',project_name)
  )
  
}
pc_genes <- readLines('pc_genes.txt')

#Load the main object

outdir_obj <- 'Analysis/Objects/ProteinCoding/'
outdir_app <- 'Analysis/Apps/ProteinCoding/'
so <- readRDS('LottieGreyWilson-RAW-noMTRP.rds')
counts <- GetAssayData(so, assay = "RNA")
counts <- counts[which(rownames(counts) %in% pc_genes),]
so <- subset(so, features = rownames(counts))
ncores <- 105
library_to_id_map <- c('C1' = 'D6',
                       'D1' = 'D10',
                       'E1' = 'D18',
                       'F1' = 'D23',
                       'G8' = 'P0',
                       'H8' = 'P2',
                       'A2' = 'LBO',
                       'A1' = 'Carola_14w',
                       'B1' = 'Carola_12w',
                       'G1' = 'Carola_5w',
                       'H1' = 'Carola_20w',
                       'B8' = 'iHB0_5-7',
                       'C8' = 'iHB0_Hc',
                       'D8' = 'iHB0_Chl')
so@meta.data$ID <- library_to_id_map[so@meta.data$library]

doublets1 <- read.csv('doublets_doubletdetection.csv')
doublets1 <- doublets1 %>%
  mutate(Doublet_dd = if_else(Doublet_dd == 0, 'False', 'True'))
doublets1$Barcode <- sapply(strsplit(doublets1$bc, "_"), function(x) x[2])
doublets1$bc <- NULL

doublets2 <- read.csv('doublets_scrublet.csv')
doublets2$Barcode <- sapply(strsplit(doublets2$bc, "_"), function(x) x[2])
doublets2$bc <- NULL
doublets <- merge(doublets1,doublets2,by='Barcode')
doublets$Doublet <- ifelse(doublets$Doublet_dd == "True" | doublets$Doublet_scrublet == "True", "yes", "no")
rownames(doublets) <- doublets$Barcode
doublets <- doublets['Doublet']

rows <- rownames(so@meta.data)
so@meta.data <- merge(so@meta.data, doublets, by = 0)
rownames(so@meta.data) <- rows
so@meta.data$Row.names <- NULL


D6 <- subset(x = so, idents = 'C1')
makeCA(D6,ncores,outdir_obj,outdir_app,'SCT','D6')
rm(D6)
gc()

D10 <- subset(x = so, idents = 'D1')
makeCA(D10,ncores,outdir_obj,outdir_app,'SCT','D10')
rm(D10)
gc()

D18 <- subset(x = so, idents = 'E1')
makeCA(D18,ncores,outdir_obj,outdir_app,'SCT','D18')
rm(D18)
gc()

D23 <- subset(x = so, idents = 'F1')
makeCA(D23,ncores,outdir_obj,outdir_app,'SCT','D23')
rm(D23)
gc()

P0 <- subset(x = so, idents = 'G8')
makeCA(P0,ncores,outdir_obj,outdir_app,'SCT','P0')
rm(P0)
gc()

P2 <- subset(x = so, idents = 'H8')
makeCA(P2,ncores,outdir_obj,outdir_app,'SCT','P2')
rm(P2)
gc()

LBO <- subset(x = so, idents = 'A2')
makeCA(LBO,ncores,outdir_obj,outdir_app,'SCT','LBO')
rm(LBO)
gc()

twoD <- subset(x = so, idents = c('C1','D1','E1','F1'))
makeCA(twoD,ncores,outdir_obj,outdir_app,'SCT','2D-Timecourse')
rm(twoD)
gc()

threeD <- subset(x = so, idents = c('G8','H8','A2'))
makeCA(threeD,ncores,outdir_obj,outdir_app,'SCT','3D-Timecourse')
rm(threeD)
gc()

iHBO <- subset(x = so, idents = c('B8','C8','C8'))
makeCA(iHBO,ncores,outdir_obj,outdir_app,'SCT','iHBO-Timecourse')
rm(iHBO)
gc()





