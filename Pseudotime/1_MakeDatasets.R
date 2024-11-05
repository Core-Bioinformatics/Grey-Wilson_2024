options(warn=1)
library(Seurat)
library(dplyr)
library(ClustAssess)

makeCA <- function(so,ncores,outdir_obj,assay_name='SCT',project_name){
  print(project_name)
  
  var_features <- so@assays[[assay_name]]@var.features
  n_abundant <- 3000
  most_abundant_genes <- rownames(so@assays[[assay_name]])[order(Matrix::rowSums(so@assays[[assay_name]]),
                                                                 decreasing = TRUE
  )][1:n_abundant]
  
  steps = seq(from = 1000, to = mean(so@meta.data$nFeature_RNA), by = 250)
  steps = steps[-length(steps)]
  
  gc()
  
  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )

  doParallel::registerDoParallel(cl = my_cluster)
  
  RhpcBLASctl::blas_set_num_threads(ncores)
  clustassess_autom <- automatic_stability_assessment(
    expression_matrix = as.matrix(so@assays$SCT$scale.data),
    n_repetitions = 8,
    n_neigh_sequence = seq(from = 5, to = 50, by = 5),
    resolution_sequence = seq(from = 0.1, to = 1, by = 0.1),
    features_sets = list(
      "HV" = var_features,
      "MA" = most_abundant_genes
    ),
    steps = list(
      "HV" = steps,
      "MA" = steps
    ),
    n_top_configs = 2,
    umap_arguments = list(
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine"
    ),
    save_temp = FALSE,
    verbose = TRUE
  )
  
  
  
  saveRDS(clustassess_autom, paste0(outdir_obj,project_name,'-2.rds')) 
  parallel::stopCluster(cl = my_cluster)
  
  ClustAssess::write_shiny_app(
    object = so,
    assay_name = "SCT",
    clustassess_object = clustassess_autom,
    project_folder = paste0(outdir_obj,project_name,'-2/'),
    shiny_app_title = paste0('Lottie - ',project_name)
  )
  
}

outdir_obj <- 'ca/'
ncores <- 8
options(warn=-1)


#now 2d of selected clusters in lottes ppt
so <- readRDS('2D-Timecourse-so-annotated-removeCells.rds')
Idents(so) <- so@meta.data$Annotations
so <- subset(so, idents = c('Early HB','Early HPB','Hepatic endoderm',
                            'Hepatic Endonderm','Posterior foregut',
                            'Proliferating posterior foregut',
                            'Proliferative cells to split'), invert = FALSE)
gc()
so <- SCTransform(so, vst.flavor = "v2", return.only.var.genes=FALSE, verbose=TRUE)
saveRDS(so,'2D_clusterSubset.rds')

so <- readRDS('2D_clusterSubset.rds')
makeCA(so,ncores,outdir_obj,'SCT','2D_clusterSubset')
