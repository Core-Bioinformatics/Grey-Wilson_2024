library(Seurat)
library(gprofiler2)
library(ggplot2)
library(reshape2)
library(viridis)
so <- readRDS('/Volumes/T7Shield/Lottie/Analysis/ProteinCoding/Annotated/May2024Annotations/2D-Timecourse.rds')

Idents(so) <- so@meta.data$Annotations
CMMarkers <- FindMarkers(so,ident.1 = 'Cardiomyocytes', logfc.threshold=0.5)
sigCM <- rownames(CMMarkers)[CMMarkers$p_val_adj<0.05]
EndMarkers <- FindMarkers(so,ident.1 = 'Endothelial', logfc.threshold=0.5)
sigEnd <- rownames(EndMarkers)[EndMarkers$p_val_adj<0.05]

allGenes <- rownames(so)


resultsCM <- gost(query = sigCM, 
                organism = "hsapiens", 
                significant = TRUE, 
                user_threshold = 0.05, 
                correction_method = "fdr",
                custom_bg = allGenes)$result

resultsCM <- apply(resultsCM,2,as.character)
write.table(resultsCM, file = '/Volumes/T7Shield/Lottie/Analysis/MarkersPerTimePoint/2D-CM_GEO.tsv', quote = FALSE, row.names = FALSE, sep = '\t')

resultsEnd <- gost(query = sigEnd, 
                   organism = "hsapiens", 
                   significant = TRUE, 
                   user_threshold = 0.05, 
                   correction_method = "fdr",
                   custom_bg = allGenes)$result

resultsEnd <- apply(resultsEnd,2,as.character)
write.table(resultsEnd, file = '/Volumes/T7Shield/Lottie/Analysis/MarkersPerTimePoint/2D-End_GEO.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
resultsCM <- as.data.frame(resultsCM)
resultsEnd <- as.data.frame(resultsEnd)

resultsCMfilter <- resultsCM[resultsCM$term_id %in% c('GO:0009653','GO:0048731','GO:0035295','GO:0072359','GO:0001944','GO:0001568','GO:0001525','GO:0030154','GO:0007507','GO:0003013'),]
resultsEndfilter <- resultsEnd[resultsEnd$term_id %in% c('GO:0009653','GO:0072359','GO:0048731','GO:0030029','GO:0030016','GO:0043292','GO:0030017','GO:0007507','GO:0060047','GO:0061061','GO:0003012','GO:0003205','GO:0014706','GO:0048738'),]


resultsCMfilter$ratio <- as.numeric(resultsCMfilter$intersection_size) / as.numeric(resultsCMfilter$term_size)
resultsCMfilter$p_value <- as.numeric(resultsCMfilter$p_value)


p <- ggplot(resultsCMfilter, aes(x = ratio, y = reorder(term_name, ratio))) +
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Enriched terms Cardiomyocytes",
    x = "Ratio",
    y = "Terms",
    color = "-log10(Pvalue)",
    size = "Number"
  ) +
  guides(size = "none") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
ggsave("/Volumes/T7Shield/Lottie/Analysis/Plots/GEO_CM_basic.svg", plot = p, device = "svg", width = 10, height = 10)

resultsEndfilter$ratio <- as.numeric(resultsEndfilter$intersection_size) / as.numeric(resultsEndfilter$term_size)
resultsEndfilter$p_value <- as.numeric(resultsEndfilter$p_value)


p <- ggplot(resultsEndfilter, aes(x = ratio, y = reorder(term_name, ratio))) +
  geom_point(aes(size = intersection_size, color = -log10(p_value))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Enriched terms Endothelial",
    x = "Ratio",
    y = "Terms",
    color = "-log10(Pvalue)",
    size = "Number"
  ) +
  guides(size = "none") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
ggsave("/Volumes/T7Shield/Lottie/Analysis/Plots/GEO_EN_basic.svg", plot = p, device = "svg", width = 10, height = 10)
#####################
#This analisys is not wanted
comp <- 'Cardiomyocytes'
Idents(so) <- so@meta.data$Annotations
rest <- unique(so@meta.data$Annotations)
rest <-  rest[rest != comp]
CMterms <- c('GO:0009653','GO:0048731','GO:0035295','GO:0072359','GO:0001944','GO:0001568','GO:0001525','GO:0030154','GO:0007507','GO:0003013')
CMterms_names <- c('Anatomical structure morphogenesis (GO:0009653)',
                   'System development (GO:0048731)',
                   'Tube development (GO:0035295)',
                   'Circulatory system development (GO:0072359)',
                   'Vasculature development (GO:0001944)',
                   'Blood vessel development (GO:0001568)',
                   'Angiogenesis (GO:0001525)',
                   'Cell differentiation (GO:0030154)',
                   'Heart development (GO:0007507)',
                   'Circulatory system process (GO:0003013)')

CMterms <- c('GO:0072359','GO:0030016','GO:0043292','GO:0030017','GO:0007507','GO:0061061','GO:0003205','GO:0014706','GO:0006936','GO:0031674','GO:0030018','GO:0042692')
CMterms_names <- c('Circulatory system development (GO:0072359)',
                   'Myofibril (GO:0030016)',
                   'Contractile fiber (GO:0043292)',
                   'Sarcomere (GO:0030017)',
                   'Heart development (GO:0007507)',
                   'Muscle structure development (GO:0061061)',
                   'Cardiac chamber development (GO:0003205)',
                   'Striated muscle tissue development (GO:0014706)',
                   'Muscle contraction (GO:0006936)',
                   'I band (GO:0031674)',
                   'Z disc (GO:0030018)',
                   'Muscle cell differentiation (GO:0042692)')

allGenes <- rownames(so)

mat_p <- matrix(nrow = length(rest), ncol = length(CMterms))
colnames(mat_p) <- CMterms
rownames(mat_p) <- rest
mat_ratio <- mat_p

for (cellType in rest){
  print(cellType)
  CMMarkers <- FindMarkers(so,ident.1 = comp,ident.2 = cellType, logfc.threshold=0.5)
  sigCM <- rownames(CMMarkers)[CMMarkers$p_val_adj<0.05]
  resultsCM <- gost(query = sigCM, 
                    organism = "hsapiens", 
                    significant = TRUE, 
                    user_threshold = 0.05, 
                    correction_method = "fdr",
                    custom_bg = allGenes)$result
  
  resultsCMfilter <- resultsCM[resultsCM$term_id %in% CMterms,]
  
  missing_terms <- setdiff(CMterms, resultsCM$term_id)
  missing_df <- data.frame(term_id = missing_terms)
  
  if (length(missing_terms) > 0){
    for(col in colnames(resultsCM)[-1]) {
      missing_df[[col]] <- 0
    }
    missing_df$term_id <- missing_terms
    missing_df$query <- 0
    missing_df$p_value <- NaN
    resultsCMfilter <- rbind(resultsCMfilter, missing_df)
  }
  
  resultsCMfilter$ratio <- resultsCMfilter$intersection_size/resultsCMfilter$term_size
  resultsCMfilter <- resultsCMfilter[c('term_id','term_name','p_value','ratio')]
  p_values <- resultsCMfilter[c('term_id','p_value')]
  rownames(p_values) <- p_values$term_id
  p_values$term_id <- NULL
  p_values <- as.data.frame(t(p_values))
  p_values <- as.matrix(p_values %>% dplyr::select(all_of(CMterms)))
  mat_p[cellType,] <- p_values[1,]
  
  
  ratio <- resultsCMfilter[c('term_id','ratio')]
  rownames(ratio) <- ratio$term_id
  ratio$term_id <- NULL
  ratio <- as.data.frame(t(ratio))
  ratio <- as.matrix(ratio %>% dplyr::select(all_of(CMterms)))
  mat_ratio[cellType,] <- ratio[1,]
}


df_ratio <- as.data.frame(t(mat_ratio))
rownames(df_ratio) <- CMterms_names
df_p <- as.data.frame(t(mat_p))
rownames(df_p) <- CMterms_names

df_ratio$GO_term <- rownames(df_ratio)
df_p$GO_term <- rownames(df_p)

df_ratio_long <- melt(df_ratio, id.vars = "GO_term", variable.name = "Cell_type", value.name = "Ratio")
df_p_long <- melt(df_p, id.vars = "GO_term", variable.name = "Cell_type", value.name = "P_value")

df_combined <- merge(df_ratio_long, df_p_long, by = c("Cell_type", "GO_term"))

# Create the bubble plot
p <- ggplot(df_combined, aes(x = Cell_type, y = GO_term, size = Ratio, color = -log10(P_value))) +
  geom_point() +
  scale_color_viridis_c() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Enriched GEO terms in Cardiomyocytes",
       x = "Cell Type",
       y = "GO Term",
       size = "Ratio",
       color = "-log10(P-value)")+
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("/Volumes/T7Shield/Lottie/Analysis/Plots/GEO_CM_2_2.svg", plot = p, device = "svg", width = 10, height = 10)
  
  


comp <- 'Endothelial'
rest <- unique(so@meta.data$Annotations)
rest <-  rest[rest != comp]
ENterms <- c('GO:0009653','GO:0072359','GO:0048731','GO:0030029','GO:0030016','GO:0043292','GO:0030017','GO:0007507','GO:0060047','GO:0061061','GO:0003012','GO:0003205','GO:0014706','GO:0048738')
ENterms_names <- c('Anatomical structure morphogenesis (GO:0009653)',
                   'Circulatory system development (GO:0072359)',
                   'System development (GO:0048731)',
                   'Actin filament-based process (GO:0030029)',
                   'Myofibril (GO:0030016)',
                   'Contractile fiber (GO:0043292)',
                   'Sarcomere (GO:0030017)',
                   'Heart development (GO:0007507)',
                   'Heart contraction (GO:0060047)',
                   'Muscle structure development (GO:0061061)',
                   'Muscle system process (GO:0003012)',
                   'Cardiac chamber development (GO:0003205)',
                   'Striated muscle tissue development (GO:0014706)',
                   'Cardiac muscle tissue development (GO:0048738)')

ENterms <- c('GO:0072359','GO:0001944','GO:0001568','GO:0001525','GO:0007507')
ENterms_names <- c('Circulatory system development (GO:0072359)',
                   'Vasculature development (GO:0001944)',
                   'Blood vessel development (GO:0001568)',
                   'Angiogenesis (GO:0001525)',
                   'Heart development (GO:0007507)')
allGenes <- rownames(so)

mat_p <- matrix(nrow = length(rest), ncol = length(ENterms))
colnames(mat_p) <- ENterms
rownames(mat_p) <- rest
mat_ratio <- mat_p

for (cellType in rest){
  print(cellType)
  ENMarkers <- FindMarkers(so,ident.1 = comp,ident.2 = cellType, logfc.threshold=0.5)
  sigEN <- rownames(ENMarkers)[ENMarkers$p_val_adj<0.05]
  resultsEN <- gost(query = sigEN, 
                    organism = "hsapiens", 
                    significant = TRUE, 
                    user_threshold = 0.05, 
                    correction_method = "fdr",
                    custom_bg = allGenes)$result
  
  resultsENfilter <- resultsEN[resultsEN$term_id %in% ENterms,]
  missing_terms <- setdiff(ENterms, resultsEN$term_id)
  missing_df <- data.frame(term_id = missing_terms)
  
  if (length(missing_terms) > 0){
    for(col in colnames(resultsEN)[-1]) {
      missing_df[[col]] <- 0
    }
    missing_df$term_id <- missing_terms
    missing_df$query <- 0
    missing_df$p_value <- NaN
    resultsENfilter <- rbind(resultsENfilter, missing_df)
  }
  resultsENfilter$ratio <- resultsENfilter$intersection_size/resultsENfilter$term_size
  resultsENfilter <- resultsENfilter[c('term_id','term_name','p_value','ratio')]
  p_values <- resultsENfilter[c('term_id','p_value')]
  rownames(p_values) <- p_values$term_id
  p_values$term_id <- NULL
  p_values <- as.data.frame(t(p_values))
  p_values <- as.matrix(p_values %>% dplyr::select(all_of(ENterms)))
  mat_p[cellType,] <- p_values[1,]
  
  
  ratio <- resultsENfilter[c('term_id','ratio')]
  rownames(ratio) <- ratio$term_id
  ratio$term_id <- NULL
  ratio <- as.data.frame(t(ratio))
  ratio <- as.matrix(ratio %>% dplyr::select(all_of(ENterms)))
  mat_ratio[cellType,] <- ratio[1,]
}


df_ratio <- as.data.frame(t(mat_ratio))
rownames(df_ratio) <- ENterms_names
df_p <- as.data.frame(t(mat_p))
rownames(df_p) <- ENterms_names

df_ratio$GO_term <- rownames(df_ratio)
df_p$GO_term <- rownames(df_p)

df_ratio_long <- melt(df_ratio, id.vars = "GO_term", variable.name = "Cell_type", value.name = "Ratio")
df_p_long <- melt(df_p, id.vars = "GO_term", variable.name = "Cell_type", value.name = "P_value")

df_combined <- merge(df_ratio_long, df_p_long, by = c("Cell_type", "GO_term"))


p <- ggplot(df_combined, aes(x = Cell_type, y = GO_term, size = Ratio, color = -log10(P_value))) +
  geom_point() +
  scale_color_viridis_c() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Enriched GEO terms in Endothelial",
       x = "Cell Type",
       y = "GO Term",
       size = "Ratio",
       color = "-log10(P-value)") +
  theme(
    panel.background = element_blank(),      
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("/Volumes/T7Shield/Lottie/Analysis/Plots/GEO_End_2_2.svg", plot = p, device = "svg", width = 10, height = 10)
