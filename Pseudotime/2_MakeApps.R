#make  pseudotime apps
library(Seurat)
library(monocle3)
library(SeuratWrappers)
source("pseudotime_roots_generate_shiny.R")
source("pseudotime_subset_generate_shiny.R")

write_object.data_frame(
  ids_df = read.csv("pseudotime_sample_combinations5.csv", comment.char = "#"),
  seurat_ca_corresp = read.csv("pseudotime_so_CA_obj5.csv", comment.char = "#"),
  target_app_dir = "",
  use_closed_loops = FALSE,
  prefix = "",
  learn_graph_controls = list(
    eps = 1e-5,
    maxiter = 100
  ),
  nodes_per_log10_cells = 30
)


write_object(
  mon_obj = readRDS(file.path('', "monocle_object.rds")),
  trajectory_id = paste('3D', 'HB', sep = "-"),
  start_genes = c('AFP'),
  end_genes = c('APOM'),
  start_expression_thresh = 4,
  end_expression_thresh = 1.9,
  start_relax_ngenes = 0,
  end_relax_ngenes = 0,
  ncores = 2,
  learn_graph_controls = list(
    eps = 1e-5,
    maxiter = 100
  ),
  use_closed_loops = FALSE,
  use_partitions = FALSE,
  nodes_per_log10_cells = 30,
  output_dir = file.path('', '3D-HB')
)


write_object(
  mon_obj = readRDS(file.path('', "monocle_object.rds")),
  trajectory_id = paste('2D', 'Subset', sep = "-"),
  start_genes = c('FOXA2'),
  end_genes = c('AFP'),
  start_expression_thresh = 1.2,
  end_expression_thresh = 4.5,
  start_relax_ngenes = 0,
  end_relax_ngenes = 0,
  ncores = 2,
  learn_graph_controls = list(
    eps = 1e-5,
    maxiter = 100
  ),
  use_closed_loops = FALSE,
  use_partitions = FALSE,
  nodes_per_log10_cells = 30,
  output_dir = file.path('', '2D-Subset')
)