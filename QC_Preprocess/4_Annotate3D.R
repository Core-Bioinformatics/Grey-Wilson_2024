library(dplyr)
library(Seurat)
library(ShinyCell)
library(ClustAssess)
####
so <- readRDS('3D-Timecourse-so.rds')
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ca <- readRDS('3D-Timecourse-ca.rds')
clusters_21 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`21`$partitions[[1]]$mb
ecc_21 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`21`$ecc

clusters_26 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`26`$partitions[[1]]$mb
ecc_26 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`26`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_21 <- ecc_21
embedding$clusters_21_SLM <- as.factor(clusters_21)
embedding$ecc_26 <- ecc_26
embedding$clusters_26_SLM <- as.factor(clusters_26)

embedding <- embedding[c('ecc_21','clusters_21_SLM','ecc_26','clusters_26_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)

so@meta.data$Annotations[so@meta.data$clusters_21_SLM==1] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==2] <- 'late HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==3] <- 'late HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==4] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==5] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==6] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==7] <- 'late HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==8] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==9] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==10] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==11] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==12] <- 'GBC'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==13] <- 'Early HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==14] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==15] <- 'Bridge pop. (diff. HPB-to-HB)'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==16] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==17] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==18] <- 'late HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==19] <- 'HPB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==20] <- 'late HB'
so@meta.data$Annotations[so@meta.data$clusters_21_SLM==21] <- '21 - query Cells'

#Now Run CA again

project_name <- '3D-Timecourse'
Idents(so) <- so@meta.data$clusters_21_SLM
#Remove 21
so <- subset(so,idents=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
so <- SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)
saveRDS(so,'3D-so-FinalNoClustering.rds')

clustassess_autom <-readRDS('3D-ca-FinalNoClustering.rds')
#18 clusters
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
clusters_18 <- clustassess_autom$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`18`$partitions[[1]]$mb
ecc_18 <- clustassess_autom$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`18`$ecc

embedding <- data.frame(clustassess_autom$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_18 <- ecc_18
embedding$clusters_18_SLM <- as.factor(clusters_18)

embedding <- embedding[c('ecc_18','clusters_18_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)

scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "3D-Timecourse",
                        shiny.title='3D-Timecourse')


#Now I need to do the opposite, from the 3D, annotate LBO,P0,P2
so <- readRDS('3D-so-FinalNoClustering.rds')
anno <- so@meta.data[c('Annotations')]
so <-  readRDS('LBO-so.rds')
meta <- so@meta.data
meta <- merge(meta,anno,by='row.names', all.x = TRUE)
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL
so@meta.data <- meta

ca <- readRDS('LBO-ca.rds')
clusters_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$partitions[[1]]$mb
ecc_7 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`7`$ecc

clusters_11 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`11`$partitions[[1]]$mb
ecc_11 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`11`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_7 <- ecc_7
embedding$clusters_7_SLM <- as.factor(clusters_7)
embedding$ecc_11 <- ecc_11
embedding$clusters_11_SLM <- as.factor(clusters_11)

embedding <- embedding[c('ecc_7','clusters_7_SLM','ecc_11','clusters_11_SLM')]
embedding <- embedding[rownames(embedding) %in% rownames(so@meta.data), ]
so@meta.data <- cbind(so@meta.data,embedding)
saveRDS(so,'LBO-so-annotated.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "LBO-Annotated",
                        shiny.title='LBO')
#P0
so <-  readRDS('P0-so.rds')
meta <- so@meta.data
meta <- merge(meta,anno,by='row.names', all.x = TRUE)
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL
so@meta.data <- meta

ca <- readRDS('P0-ca.rds')
clusters_5 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`5`$partitions[[1]]$mb
ecc_5 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`5`$ecc

clusters_10 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`10`$partitions[[1]]$mb
ecc_10 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`10`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_5 <- ecc_5
embedding$clusters_5_SLM <- as.factor(clusters_5)
embedding$ecc_10 <- ecc_10
embedding$clusters_10_SLM <- as.factor(clusters_10)

embedding <- embedding[c('ecc_5','clusters_5_SLM','ecc_10','clusters_10_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)
saveRDS(so,'P0-so-annotated.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "P0-Annotated",
                        shiny.title='P0')
#P2
so <-  readRDS('P2-so.rds')
meta <- so@meta.data
meta <- merge(meta,anno,by='row.names', all.x = TRUE)
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL
so@meta.data <- meta
ca <- readRDS('P2-ca.rds')
clusters_4 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`4`$partitions[[1]]$mb
ecc_4 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`4`$ecc

clusters_12 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`12`$partitions[[1]]$mb
ecc_12 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`12`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_4 <- ecc_4
embedding$clusters_4_SLM <- as.factor(clusters_4)
embedding$ecc_12 <- ecc_12
embedding$clusters_12_SLM <- as.factor(clusters_12)

embedding <- embedding[c('ecc_4','clusters_4_SLM','ecc_12','clusters_12_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)
saveRDS(so,'P2-so-annotated.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "P2-Annotated",
                        shiny.title='P2')
