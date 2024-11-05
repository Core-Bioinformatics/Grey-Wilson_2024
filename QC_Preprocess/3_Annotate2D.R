library(dplyr)
library(Seurat)
library(ShinyCell)
library(ClustAssess)
#First Annotate the individual timepoints for 2D
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#D6
so <- readRDS('Analysis/Objects/ProteinCoding/D6-so.rds')
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ca <- readRDS('Analysis/Objects/ProteinCoding/D6-ca.rds')
clusters_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$partitions[[1]]$mb
ecc_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$ecc

clusters_15 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`15`$partitions[[1]]$mb
ecc_15 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`15`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_7 <- ecc_7
embedding$clusters_7_SLM <- as.factor(clusters_7)
embedding$ecc_15 <- ecc_15
embedding$clusters_15_SLM <- as.factor(clusters_15)

embedding <- embedding[c('ecc_7','clusters_7_SLM','ecc_15','clusters_15_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==1] <- 'Posterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==2] <- 'Proliferating posterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==3] <- 'Intermediate'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==4] <- 'Posterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==5] <- 'Anterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==6] <- 'Hepatic Endonderm'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==7] <- 'Endothelial'

saveRDS(so,'D6-so-Final.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "D6",
                        shiny.title='D6')
#D10
#K6 from D10 is removed
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so <- readRDS(' D10-so.rds')
ca <- readRDS(' D10-ca.rds')
clusters_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$partitions[[1]]$mb
ecc_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_7 <- ecc_7
embedding$clusters_7_SLM <- as.factor(clusters_7)

embedding <- embedding[c('ecc_7','clusters_7_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)


so@meta.data$Annotations[so@meta.data$clusters_7_SLM==1] <- 'Hepatic endoderm'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==2] <- 'Anterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==3] <- 'Airway foregut'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==4] <- 'Endothelial'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==7] <- 'Cardiomyocytes'
so@meta.data$Annotations[so@meta.data$clusters_7_SLM==5] <- 'Proliferative cells to split'
#Remove Cluster 6
Idents(so) <- so@meta.data$clusters_7_SLM
so <- subset(so,idents=c(1,2,3,4,5,7))
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
so <- SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)
clustassess_autom <- readRDS('D10-removeK6-ca.rds')
embedding <- data.frame(clustassess_autom$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))


scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "D10-removeK6",
                        shiny.title='D10-removeK6')
#D18
so <- readRDS(' D18-so.rds')
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ca <- readRDS(' D18-ca.rds')
clusters_9 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`9`$partitions[[1]]$mb
ecc_9 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`9`$ecc

clusters_13 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`13`$partitions[[1]]$mb
ecc_13 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`13`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_9 <- ecc_9
embedding$clusters_9_SLM <- as.factor(clusters_9)
embedding$ecc_13 <- ecc_13
embedding$clusters_13_SLM <- as.factor(clusters_13)

embedding <- embedding[c('ecc_9','clusters_9_SLM','ecc_13','clusters_13_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)

so@meta.data$Annotations[so@meta.data$clusters_9_SLM==1] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==2] <- 'Anterior Foregut'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==3] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==4] <- 'Airway Foregut'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==5] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==6] <- 'Endothelial'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==7] <- 'Cardiomyocytes'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==8] <- 'Trachea Progenitor'
so@meta.data$Annotations[so@meta.data$clusters_9_SLM==9] <- 'Anterior foregut'


saveRDS(so,'D18-so-Final.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "",
                        shiny.title='D18')
#D23
so <- readRDS('D23-so.rds')
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ca <- readRDS('D23-ca.rds')
clusters_11 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`11`$partitions[[1]]$mb
ecc_11 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`11`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_11 <- ecc_11
embedding$clusters_11_SLM <- as.factor(clusters_11)

embedding <- embedding[c('ecc_11','clusters_11_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)



so@meta.data$Annotations[so@meta.data$clusters_11_SLM==1] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==2] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==3] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==4] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==5] <- 'Trachea Progenitor'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==6] <- 'Anterior foregut'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==7] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==8] <- 'Early HB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==9] <- 'Early HPB'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==10] <- 'Cardiomyocytes'
so@meta.data$Annotations[so@meta.data$clusters_11_SLM==11] <- 'Endothelial'


saveRDS(so,'D23-so-Final.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "D23",
                        shiny.title='D23')

#Now Reverse annotate the 2D overall timepoint

so <- readRDS('D6-so-Final.rds')
D6 <- so@meta.data[c('Annotations')]
so <- readRDS('D10-removeK6-so.rds')
D10 <- so@meta.data[c('Annotations')]
so <- readRDS('D18-so-Final.rds')
D18 <- so@meta.data[c('Annotations')]
so <- readRDS('D23-so-Final.rds')
D23 <- so@meta.data[c('Annotations')]

anno <- rbind(D6,D10,D18,D23)
so <- readRDS('2D-Timecourse-so.rds')
meta <- so@meta.data
meta <- merge(meta,anno,by = "row.names", all = TRUE)
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL
so@meta.data <- meta
saveRDS(so,'2D-Timecourse-so-annotated.rds')
so <- readRDS('2D-Timecourse-so-annotated.rds')
ca <- readRDS('2D-Timecourse-ca.rds')
clusters_8 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`8`$partitions[[1]]$mb
ecc_8 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`8`$ecc

clusters_20 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`26`$partitions[[1]]$mb
ecc_20 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`26`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_8 <- ecc_8
embedding$clusters_8_SLM <- as.factor(clusters_8)
embedding$ecc_20 <- ecc_20
embedding$clusters_20_SLM <- as.factor(clusters_20)

embedding <- embedding[c('ecc_8','clusters_8_SLM','ecc_20','clusters_20_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)
Idents(so) <- so@meta.data$Annotations
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = '2D-Annotated',
                        shiny.title='2D-Annotated')


