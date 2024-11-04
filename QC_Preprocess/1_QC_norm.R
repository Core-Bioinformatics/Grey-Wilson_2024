library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
library(readr)
library(tibble)
library(ClustAssess)
library(gridExtra)
set.seed(42)

setwd('')

create.so = function(filenames, super.folder='counts', min.cells=0, min.features=0, project.name='LottieGreyWilson'){
  mat = list()
  lib.name = NULL
  for (i in 1:length(folders)){
    current.file = filenames[i]
    #print(paste0(super.folder, '/', current.file, '/filtered_feature_bc_matrix'))
    temp.mat = Read10X(paste0(super.folder, '/', current.file))
    mat[[i]] = temp.mat
    filename_stripped = strsplit(current.file, split='.h5', fixed=TRUE)[[1]][1]
    lib.name = c(lib.name, rep(filename_stripped, ncol(temp.mat)))
  }
  cts = do.call(cbind, mat)
  cts = cts[Matrix::rowSums(cts) > 0,]
  colnames(cts) = paste0(lib.name, colnames(cts))
  
  # create metadata
  meta = data.frame(library=lib.name)
  rownames(meta) = colnames(cts)
  
  so = CreateSeuratObject(counts=cts, min.cells=min.cells, min.features=min.features, meta.data=meta, project=project.name)
  return(so)
}

min.cells = 0
min.features = 0
folders = list.files('cr-output/')
so = create.so(folders, super.folder='cr-output/')

# filter cells
Idents(so) = so@meta.data$library
mt.genes=grep("^MT-", rownames(so), value=FALSE, ignore.case=TRUE)
rp.genes=grep("^RP[SL]", rownames(so), value=FALSE, ignore.case=TRUE)
so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[rp.genes])

pdf('Analysis/Plots/soQC.pdf')
ggplot(so@meta.data, aes(x=library)) + geom_bar() + theme(axis.text.x=element_text(angle=45, hjust=1))
VlnPlot(so, features = c("nFeature_RNA"), ncol = 1, pt.size=0,raster=FALSE)
VlnPlot(so, features = c("nCount_RNA"), ncol = 1, pt.size=0,raster=FALSE)
VlnPlot(so, features = c("nCount_RNA"), ncol = 1, pt.size=0,raster=FALSE,log=T)
VlnPlot(so, features = c("percent.mt"), ncol = 1, pt.size=0,raster=FALSE)
VlnPlot(so, features = c("percent.mt"), ncol = 1, pt.size=0,raster=FALSE,log=T)
VlnPlot(so, features = c("percent.rp"), ncol = 1, pt.size=0,raster=FALSE)
VlnPlot(so, features = c("percent.rp"), ncol = 1, pt.size=0,raster=FALSE,log=T)
ggplot(so@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(alpha=0.03)
ggplot(so@meta.data, aes(x=percent.mt, y=percent.rp)) + geom_point(alpha=0.03)                                                                    
ggplot(so@meta.data, aes(x=percent.mt, y=nCount_RNA)) + geom_point(alpha=0.03)
ggplot(so@meta.data, aes(x=percent.rp, y=nCount_RNA)) + geom_point(alpha=0.03)
dev.off()


#Perform filterning
length(rownames(so@meta.data))
#157839
so = subset(so, nFeature_RNA>500 & nFeature_RNA<5000 & nCount_RNA>2e3 & nCount_RNA<2e4 & percent.mt<5 & percent.rp<10 )
length(rownames(so@meta.data))
#107326
so = so[-c(mt.genes, rp.genes),]
so[['raw.seq.depth.no.MTRP']] = Matrix::colSums(GetAssayData(so, assay='RNA', slot='counts'))
saveRDS(so,'Analysis/Objects/LottieGreyWilson-RAW-noMTRP.rds')

# Separate the samples into experiment
so <- readRDS('Analysis/Objects/LottieGreyWilson-RAW-noMTRP.rds')
so.timeSeries <- subset(x = so, idents = c('C1','D1','E1','F1','G8','H8','A2'))
so.differentiation <- subset(x = so, idents = c('B8','C8','D8'))
so.Carola <- subset(x = so, idents = c('A1','B1','G1','H1'))

saveRDS(so.timeSeries,'Analysis/Objects/timeSeries-RAW-noMTRP.rds')
saveRDS(so.differentiation,'Analysis/Objects/differentiation-RAW-noMTRP.rds')
saveRDS(so.Carola,'Analysis/Objects/Carola-RAW-noMTRP.rds')


obj_names <- c('timeSeries','differentiation','Carola')
objects <- c(so.timeSeries,so.differentiation,so.Carola)


so.timeSeries <-SCTransform(so.timeSeries, return.only.var.genes=FALSE, verbose=FALSE)
so.differentiation <-SCTransform(so.differentiation, return.only.var.genes=FALSE, verbose=FALSE)
so.Carola <-SCTransform(so.Carola, return.only.var.genes=FALSE, verbose=FALSE)

#I want to add some metadata
library_to_id_map <- c('C1' = 'D6',
                       'D1' = 'D10',
                       'E1' = 'D18',
                       'F1' = 'D23',
                       'G8' = 'P0',
                       'H8' = 'P2',
                       'A2' = 'LBO')
so.timeSeries@meta.data$ID <- library_to_id_map[so.timeSeries@meta.data$library]

saveRDS(so.timeSeries,'Analysis/Objects/timeSeries-SCTransformed-noMTRP.rds')
saveRDS(so.differentiation,'Analysis/Objects/differentiation-SCTransformed-noMTRP.rds')
saveRDS(so.Carola,'Analysis/Objects/Carola-SCTransformed-noMTRP.rds')


  
  