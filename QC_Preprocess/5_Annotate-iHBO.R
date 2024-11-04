so <- readRDS('Analysis/ProteinCoding/ProteinCoding/iHBO-Timecourse-so.rds')
ca <- readRDS('Lottie/Analysis/ProteinCoding/ProteinCoding/iHBO-Timecourse-ca.rds')
clusters_4 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`4`$partitions[[1]]$mb
ecc_4 <- ca$Most_Abundant$`1750`$clustering_stability$split_by_k$SLM$`4`$ecc

embedding <- data.frame(ca$Most_Abundant$`1750`$umap)
so@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(so))

embedding$ecc_4 <- ecc_4
embedding$clusters_4_SLM <- as.factor(clusters_4)


embedding <- embedding[c('ecc_4','clusters_4_SLM')]
so@meta.data <- cbind(so@meta.data,embedding)
saveRDS(so,'Analysis/ProteinCoding/Annotated/May2024Annotations/iHBO-Timecourse.rds')
scConf = ShinyCell::createConfig(so)
ShinyCell::makeShinyApp(so, scConf, gene.mapping = TRUE, gex.assay = "SCT",
                        shiny.dir = "Analysis/Apps/ProteinCoding/AnnotatedShinyCell/FinalAnnotationMay/iHBO-Timecourse",
                        shiny.title='iHBO-Timecourse')