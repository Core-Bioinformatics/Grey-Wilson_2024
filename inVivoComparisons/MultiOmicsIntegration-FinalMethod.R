library(Seurat)
library(rhdf5)
library(dplyr)
library(ggplot2)

out.dir <- 'MutiOmicsIntegration/'

so <- readRDS('so_merged_100k_without_14w_new_MA-clustassess.rds')
#I think the annotations are on this app so_merged_100k_without_14w_new_MA

meta <- as.data.frame(readRDS('so_merged_100k_without_14w_new_MA/sc1meta.rds'))
meta <- meta[c('sampleID','celltypes','UMAP_1','UMAP_2')]
rownames(meta) <- meta$sampleID

someta <- so@meta.data
someta <- merge(someta,meta,by= 'row.names')
rownames(someta) <- someta$Row.names
someta$Row.names <- NULL
someta$sampleID <- NULL
so@meta.data <- someta
DimPlot(so,group.by = 'celltypes')

pc_genes <- readLines('pcGenes.txt')
counts <- GetAssayData(so, assay = "RNA")
counts <- counts[which(rownames(counts) %in% pc_genes),]
so <- subset(so, features = rownames(counts))
cells <-rownames(so@meta.data[so@meta.data$celltypes=='hepatocytes',])
# for some reason I cant subset cells so
exp <- as.matrix(so@assays$RNA$counts)
expHep <- exp[,cells]
someta$bc <- rownames(someta)
metaHep <- someta[someta$bc%in%cells,]
soHep <- CreateSeuratObject(expHep,)
soHep@meta.data$age <- metaHep$age
soHep <- SCTransform(soHep, return.only.var.genes=FALSE, verbose=T)
exp <- as.matrix(soHep@assays$SCT@scale.data)

meta <-soHep@meta.data[c('age')]
meta$cell <- rownames(meta)
top_percent <- function(mat, rownames) {
  valid_rownames <- rownames %in% rownames(mat)
  if (any(!valid_rownames)) {
    print("Some gene names are not present in ref.")
  }
  
  filtered_rownames <- rownames[valid_rownames]
  
  result <- matrix(FALSE, nrow = length(filtered_rownames), ncol = ncol(mat))
  rownames(result) <- filtered_rownames
  colnames(result) <- colnames(mat)
  
  for (i in seq_along(filtered_rownames)) {
    row <- mat[filtered_rownames[i], ]
    threshold <- quantile(row, 0.75)  
    result[i, ] <- row >= threshold  
  }
  
  return(result)
}
bottom_percent <- function(mat, rownames) {
  valid_rownames <- rownames %in% rownames(mat)
  if (any(!valid_rownames)) {
    print("Some gene names are not present in ref.")
  }
  
  filtered_rownames <- rownames[valid_rownames]
  
  result <- matrix(FALSE, nrow = length(filtered_rownames), ncol = ncol(mat))
  rownames(result) <- filtered_rownames
  colnames(result) <- colnames(mat)
  
  for (i in seq_along(filtered_rownames)) {
    row <- mat[filtered_rownames[i], ]
    threshold <- quantile(row, 0.25) 
    result[i, ] <- row <= threshold  
  }
  
  return(result)
}
so <- readRDS('3D-Timecourse.rds')
Idents(so) <- so@meta.data$Annotations
markers <- FindAllMarkers(so, logfc.threshold=0.1, min.pct=0.4, verbose=T)

df <- data.frame(row.names = names(table(meta$age)))
timepoints <- c('late HB','Early HB')
for (t in timepoints){
  print(t)
  
  markers_sub <- markers[markers$cluster == t & markers$p_val_adj<0.05,]
  genes_total <- length(unique(markers_sub$gene))
  genes_up <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC > 0]
  genes_down <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC < 0]
  
  up <- top_percent(exp, genes_up)
  down <- bottom_percent(exp, genes_down)
  
  cell_marker_count <- (colSums(up)+colSums(down))
  selected_cells <- names(cell_marker_count[cell_marker_count > (genes_total*0.3)])
  
  selected_cell_ages <- meta %>% filter(cell %in% selected_cells)
  df[,t] <- table(selected_cell_ages$age)
  df$Total <- table(soHep@meta.data$age)
  df_proportions <- df/df$Total
  df_proportions <- df_proportions %>% select(-contains("Var1"))
  df_proportions$Total.Freq <- NULL
  df_proportions <- df_proportions[!(row.names(df_proportions) %in% c('14w')),]
  df_proportions$age <-rownames(df_proportions)
  df_mat <- reshape2::melt(df_proportions)
}

p <- ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(variable), as.factor(age))) +
        ggplot2::geom_tile(ggplot2::aes(fill = value)) +
        ggplot2::geom_text(ggplot2::aes(label = round(value, 2))) +
        ggplot2::scale_fill_gradient2(
          low = scales::muted("darkred"),
          mid = "white",
          high = scales::muted("green"),
          midpoint = 0
        ) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white"),
          axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
          axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
          axis.title = ggplot2::element_text(size = 14, face = "bold"),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
          axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
        ) +
        ggplot2::xlab('') +
        ggplot2::ylab('') +
        ggplot2::labs(fill = 'Proportion of cells')+
        ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave("3DvsInVivoHep,svg", plot = p, device = "svg", width = 6, height = 9)

#Save contribuiting markers
markers_sig <- markers[markers$p_val_adj<0.05,]
write.csv(markers_sig,paste0(out.dir,'3D_HepMarkers.csv'),quote = F,row.names = F)


calculateMarkerSim <- function(expRef,so,topPercent,p,timepoints,logfc,min.pct){
  
  markers <- FindAllMarkers(so, logfc.threshold=logfc, min.pct=min.pct, verbose=T)
  
  df <- data.frame(row.names = names(table(meta$age)))
  for (t in timepoints){
    print(t)
    
    markers_sub <- markers[markers$cluster == t & markers$p_val_adj<p,]
    genes_total <- length(unique(markers_sub$gene))
    genes_up <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC > 0]
    genes_down <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC < 0]
    
    up <- top_percent(expRef, genes_up)
    down <- bottom_percent(expRef, genes_down)
    
    cell_marker_count <- (colSums(up)+colSums(down))
    selected_cells <- names(cell_marker_count[cell_marker_count > (genes_total*topPercent)])
    
    selected_cell_ages <- meta %>% filter(cell %in% selected_cells)
    df[,t] <- table(selected_cell_ages$age)
    df$Total <- table(soHep@meta.data$age)
    df_proportions <- df/df$Total
    df_proportions <- df_proportions %>% select(-contains("Var1"))
    df_proportions$Total.Freq <- NULL
    df_proportions <- df_proportions[!(row.names(df_proportions) %in% c('14w')),]
    df_proportions$age <-rownames(df_proportions)
    df_mat <- reshape2::melt(df_proportions)
  }
  return(ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(variable), as.factor(age))) +
           ggplot2::geom_tile(ggplot2::aes(fill = value)) +
           ggplot2::geom_text(ggplot2::aes(label = round(value, 2))) +
           ggplot2::scale_fill_gradient2(
             low = scales::muted("darkred"),
             mid = "white",
             high = scales::muted("green"),
             midpoint = 0
           ) +
           ggplot2::theme(
             panel.background = ggplot2::element_rect(fill = "white"),
             axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
             axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
             axis.title = ggplot2::element_text(size = 14, face = "bold"),
             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
             axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
           ) +
           ggplot2::xlab('') +
           ggplot2::ylab('') +
           ggplot2::labs(fill = 'Proportion of cells')+
           ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))
}
so <- readRDS('iHBO-Timecourse.rds')
Idents(so) <- so@meta.data$ID
timepoints <-c('iHB0_5-7','iHB0_Hc')
p <- calculateMarkerSim(exp,so,0.3,0.1,timepoints,0.5,0.5)
p
pdf(paste0(out.dir,'iHBO_InVivoHep.pdf'),height=10,width=5)
p
dev.off()
#2D vs multiomics
so <- readRDS('2D-Timecourse.rds')
Idents(so) <- so@meta.data$Annotations
timepoints <-c('Posterior foregut','Proliferative posterior foregut','Hepatic endoderm (D6)',
               'Hepatic endoderm (D10)','Early HPB','Early HB')
#Remember to change to the top 20% instead of 25
p <- calculateMarkerSim(exp,so,0.2,0.001,timepoints,0.2,0.7)
p
pdf(paste0(out.dir,'2D_InVivoHep.pdf'),height=10,width=5)
p
dev.off()

#iHBO to organoids
so <- readRDS('iHBO-Timecourse.rds')
Idents(so) <- so@meta.data$ID
timepoints <-c('iHB0_5-7','iHB0_Hc','iHB0_Chl')
calculateMarkerSimMod <- function(expRef,so,topPercent,p,timepoints,logfc,min.pct){
  
  markers <- FindAllMarkers(so, logfc.threshold=logfc, min.pct=min.pct, verbose=T)
  
  df <- data.frame(row.names = names(table(meta$Annotations)))
  for (t in timepoints){
    print(t)
    
    markers_sub <- markers[markers$cluster == t & markers$p_val_adj<p,]
    genes_total <- length(unique(markers_sub$gene))
    genes_up <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC > 0]
    genes_down <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC < 0]
    
    up <- top_percent(expRef, genes_up)
    down <- bottom_percent(expRef, genes_down)
    
    cell_marker_count <- (colSums(up)+colSums(down))
    selected_cells <- names(cell_marker_count[cell_marker_count > (genes_total*topPercent)])
    
    selected_cell_ages <- meta %>% filter(cell %in% selected_cells)
    df[,t] <- table(selected_cell_ages$Annotations)
    df$Total <- table(InVivo@meta.data$Annotations)
    df_proportions <- df/df$Total
    df_proportions <- df_proportions %>% select(-contains("Var1"))
    df_proportions$Total.Freq <- NULL
    #df_proportions <- df_proportions[!(row.names(df_proportions) %in% c('14w')),]
    df_proportions$Annotations <-rownames(df_proportions)
    df_mat <- reshape2::melt(df_proportions)
  }
  return(ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(variable), as.factor(Annotations))) +
           ggplot2::geom_tile(ggplot2::aes(fill = value)) +
           ggplot2::geom_text(ggplot2::aes(label = round(value, 2))) +
           ggplot2::scale_fill_gradient2(
             low = scales::muted("darkred"),
             mid = "white",
             high = scales::muted("green"),
             midpoint = 0
           ) +
           ggplot2::theme(
             panel.background = ggplot2::element_rect(fill = "white"),
             axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
             axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
             axis.title = ggplot2::element_text(size = 14, face = "bold"),
             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
             axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
           ) +
           ggplot2::xlab('') +
           ggplot2::ylab('') +
           ggplot2::labs(fill = 'Proportion of cells')+
           ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)))
}
InVivo <- readRDS('InVivo-Timecourse-so-stableClusters-Annotated.rds')
expInVivoOrg <- as.matrix(InVivo@assays$RNA$counts)

meta <- InVivo@meta.data[c('Annotations')]
meta$cell <- rownames(meta)
timepoints <-c('iHB0_5-7')
p <- calculateMarkerSimMod(expInVivoOrg,so,0.4,0.05,timepoints,0.1,0.4)
ggsave("/iHBO5-7vsInVivoOrg.svg", plot = p, device = "svg", width = 6, height = 9)
timepoints <-c('iHB0_Chl')
p <- calculateMarkerSimMod(expInVivoOrg,so,0.5,0.05,timepoints,0.1,0.4)
ggsave("iHBOChlvsInVivoOrg.svg", plot = p, device = "svg", width = 6, height = 9)
timepoints <-c('iHB0_Hc')
p <- calculateMarkerSimMod(expInVivoOrg,so,0.45,0.05,timepoints,0.1,0.4)
ggsave("iHBOHcvsInVivoOrg.svg", plot = p, device = "svg", width = 6, height = 9)







topPercent_values <- seq(0.1, 0.5, by = 0.1)
p_values <- seq(0.01, 0.1, by = 0.01)
p_values <-0.05
logfc_values <- seq(0.1, 1.0, by = 0.1)
min.pct_values <- seq(0.1, 0.9, by = 0.1)



top_percent <- function(mat, rownames) {
  valid_rownames <- rownames %in% rownames(mat)
  if (any(!valid_rownames)) {
    print("Some gene names are not present in ref.")
  }
  
  filtered_rownames <- rownames[valid_rownames]
  
  result <- matrix(FALSE, nrow = length(filtered_rownames), ncol = ncol(mat))
  rownames(result) <- filtered_rownames
  colnames(result) <- colnames(mat)
  
  for (i in seq_along(filtered_rownames)) {
    row <- mat[filtered_rownames[i], ]
    threshold <- quantile(row, 0.6)  
    result[i, ] <- row >= threshold  
  }
  
  return(result)
}
bottom_percent <- function(mat, rownames) {
  valid_rownames <- rownames %in% rownames(mat)
  if (any(!valid_rownames)) {
    print("Some gene names are not present in ref.")
  }
  
  filtered_rownames <- rownames[valid_rownames]
  
  result <- matrix(FALSE, nrow = length(filtered_rownames), ncol = ncol(mat))
  rownames(result) <- filtered_rownames
  colnames(result) <- colnames(mat)
  
  for (i in seq_along(filtered_rownames)) {
    row <- mat[filtered_rownames[i], ]
    threshold <- quantile(row, 0.4) 
    result[i, ] <- row <= threshold  
  }
  
  return(result)
}
topPercent <- 0.4
topPercentOldFuns <- 0.25
pv <- 0.05
logfc <- 1
min.pct <- 0.5

calculateMarkerSim <- function(expRef,so,topPercent,p,timepoints,logfc,min.pct){
  
  markers <- FindAllMarkers(so, logfc.threshold=logfc, min.pct=min.pct, verbose=T)
  
  df <- data.frame(row.names = names(table(meta$age)))
  for (t in timepoints){
    print(t)
    
    markers_sub <- markers[markers$cluster == t & markers$p_val_adj<p,]
    genes_total <- length(unique(markers_sub$gene))
    genes_up <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC > 0]
    genes_down <- markers_sub$gene[markers_sub$cluster == t & markers_sub$avg_log2FC < 0]
    
    up <- top_percent(expRef, genes_up)
    down <- bottom_percent(expRef, genes_down)
    
    cell_marker_count <- (colSums(up)+colSums(down))
    selected_cells <- names(cell_marker_count[cell_marker_count > (genes_total*topPercent)])
    
    selected_cell_ages <- meta %>% filter(cell %in% selected_cells)
    df[,t] <- table(selected_cell_ages$age)
    df$Total <- table(soHep@meta.data$age)
    df_proportions <- df/df$Total
    df_proportions <- df_proportions %>% dplyr::select(-contains("Var1"))
    df_proportions$Total.Freq <- NULL
    df_proportions <- df_proportions[!(row.names(df_proportions) %in% c('14w')),]
    df_proportions$age <-rownames(df_proportions)
    df_mat <- reshape2::melt(df_proportions)
  }
  return(ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(variable), as.factor(age))) +
           ggplot2::geom_tile(ggplot2::aes(fill = value)) +
           ggplot2::geom_text(ggplot2::aes(label = round(value, 2))) +
           ggplot2::scale_fill_gradient2(
             low = scales::muted("darkred"),
             mid = "white",
             high = scales::muted("green"),
             midpoint = 0
           ) +
           ggplot2::theme(
             panel.background = ggplot2::element_rect(fill = "white"),
             axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
             axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
             axis.title = ggplot2::element_text(size = 14, face = "bold"),
             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
             axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
           ) +
           ggplot2::xlab('') +
           ggplot2::ylab('') +
           ggplot2::labs(fill = 'Proportion of cells')+
           ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))
}
so <- readRDS('iHBO-Timecourse.rds')
Idents(so) <- so@meta.data$ID
timepoints <-c('iHB0_5-7','iHB0_Hc')
p <- calculateMarkerSim(exp,so,topPercent,pv,timepoints,logfc,min.pct)
p
ggsave("iHBOHcvsInVivoHep.svg", plot = p, device = "svg", width = 6, height = 9)
