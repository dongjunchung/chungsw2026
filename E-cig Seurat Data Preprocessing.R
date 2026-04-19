library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(igraph)
library(readr)
library(data.table)
library(cowplot)
library(readxl)
library(rlang)
library(jsonlite)
#devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)

p1_dir <- c("./raw_feature_bc_matrix/")


sample <- Read10X( p1_dir)

save(sample, file="./beforeseurat_esig_all.RData")

load("./beforeseurat_esig_all.RData")

sample <- CreateSeuratObject(counts = sample, 
                             project = "all", 
                             min.cells = 3, #low quality genes
                             min.features = 600) #at least 200 features
length(sample$orig.ident)
length(sample1$orig.ident)
js <-  fromJSON("cells_per_tag.json", flatten=TRUE)
j <- 0
for(i in js){
  j <- j + (length(i))
}
brc <- read.csv("./Cells Per Tag.csv")
sample$patient <- "none"

sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC001"),]$Barcode] <- 1
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC002"),]$Barcode] <- 2
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC003"),]$Barcode] <- 3
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC004"),]$Barcode] <- 4
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC005"),]$Barcode] <- 5
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC006"),]$Barcode] <- 6
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC007"),]$Barcode] <- 7
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC008"),]$Barcode] <- 8
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC009"),]$Barcode] <- 17
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC010"),]$Barcode] <- 18
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC011"),]$Barcode] <- 19
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC012"),]$Barcode] <- 20
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC013"),]$Barcode] <- 21
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC014"),]$Barcode] <- 22
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC015"),]$Barcode] <- 23
sample$patient[Cells(sample) %in% brc[which(brc$Cells.Per.Tag=="BC016"),]$Barcode] <- 24

table(sample$patient)
save(sample,file="./beforeseurat_processing_esig_all.RData")
load("./beforeseurat_processing_esig_all.RData")

Idents(sample) <- "patient"
sample <- subset(sample,ident="none",invert=TRUE)
save(sample,file="./beforeseurat_processing_esig_all_remove_none.RData")
load("./beforeseurat_processing_esig_all_remove_none.RData")


#qc-check and selecting cells -filter cells that have unique feature counts over 2,500 or less than 200, We filter cells that have >5% mitochondrial counts
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

pdf("./qc_esig_all.pdf",width=40)
Idents(sample) <- "orig.ident"
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Idents(sample) <- "patient"
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1_ct <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_ct <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_ct + plot2_ct
dev.off()
#They all look fine but control and pr8 has some noise

sample <- subset(sample, subset = nFeature_RNA < 6200 & percent.mt < 5)
#integrating
save(sample, file="./esig_after_qc.Rdata")

load("./esig_after_qc.Rdata")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

# Run the standard workflow for visualization and clustering
sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
sample <- ScaleData(sample, verbose = FALSE)
sample <- RunPCA(sample, npcs = 30, verbose = FALSE)
ElbowPlot(sample)
#maybe changed from 30 to 20 to 30
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
#pca 1:30
save(sample, file ="./esig_before_findclusters.RData")
load("./esig_before_findclusters.RData")

#changed from 0.5
sample <- FindClusters(sample, resolution = 0.5)

sample$age <- 43
sample$age[which(sample$patient == 1)]<- 6
sample$age[which(sample$patient == 5)]<- 6

sample$age[which(sample$patient == 2)]<- 34
sample$age[which(sample$patient == 6)]<- 34

sample$age[which(sample$patient == 3)]<- 7
sample$age[which(sample$patient == 7)]<- 7

sample$age[which(sample$patient == 4)]<- 54
sample$age[which(sample$patient == 8)]<- 54

sample$age[which(sample$patient == 17)]<- 62
sample$age[which(sample$patient == 21)]<- 62

sample$age[which(sample$patient == 18)]<- 41
sample$age[which(sample$patient == 22)]<- 41

sample$age[which(sample$patient == 19)]<- 54
sample$age[which(sample$patient == 23)]<- 54

sample$race <- "White"
sample$race[which(sample$age == 34 | sample$age == 62)]<- "Black"
sample$race[which(sample$age == 43)] <- "Asian"

sample$sex <- "F"
sample$sex[which(sample$age == 7 | sample$age == 54  | sample$age == 43)]<- "M"

sample$disease <- "Healthy"
sample$disease[which(sample$patient == 2 | sample$patient == 6 | sample$patient == 4 | sample$patient == 8 | sample$patient == 18 | sample$patient == 22 | sample$patient == 20 | sample$patient == 24)]<- "Asthma"

sample$trt <- "Air"
sample$trt[which(sample$patient == 5 | sample$patient == 6 | sample$patient == 7 | sample$patient == 8 | sample$patient == 21 | sample$patient == 22 | sample$patient == 23 | sample$patient == 24)]<- "E-cig"


sample$sex[which(sample$sex=="F")] <- "Female"
sample$sex[which(sample$sex=="M")] <- "Male"

save(sample,file="afternorm_05_resol.RData")
load("afternorm_05_resol.RData")

counts <- sample@assays$RNA$counts
metadata <- sample@meta.data
save(counts,metadata,file = "./esig_counts.RData")

#load("./bcell_counts.RData")
# Visualization
pdf("./esig_dimplot.pdf")
p1 <- DimPlot(sample, reduction = "umap", label = TRUE, repel = TRUE)
p1

table(sample$seurat_clusters)

Idents(sample) <- "seurat_clusters"
p2 <- DimPlot(sample, reduction = "umap", split.by = "disease",label=TRUE)
p3 <- DimPlot(sample, reduction = "umap", split.by = "trt",label=TRUE)
p2/p3
dev.off()
pdf("./esig_patient.pdf",width=50)
p4 <- DimPlot(sample, reduction = "umap", split.by = "patient",label=TRUE)
p4
dev.off()

Idents(sample) <- "sex"
sample_f <- subset(sample,idents="F")
sample_m <- subset(sample,idents="M")

pdf("./various_condition_by_gender.pdf",width=10)
sample$disease_trt <- paste0(sample$disease,"_",sample$trt)
Idents(sample) <- "seurat_clusters"
p2 <- DimPlot(sample, reduction = "umap", split.by = "sex",label=TRUE)
p2

sample$sex_trt <- paste0(sample$sex,"_",sample$trt)
Idents(sample) <- "seurat_clusters"
p2 <- DimPlot(sample, reduction = "umap", split.by = "sex_trt",label=TRUE)
p2


sample$sex_dse <- paste0(sample$sex,"_",sample$disease)
Idents(sample) <- "seurat_clusters"
p3 <- DimPlot(sample, reduction = "umap", split.by = "sex_dse",label=TRUE)
p3

sample$sex_dse_trt <- paste0(sample$sex,"_",sample$disease,"_",sample$trt)

Idents(sample_m) <- "seurat_clusters"
Idents(sample_f) <- "seurat_clusters"
p3 <- DimPlot(sample_f, reduction = "umap", split.by = "sex_dse_trt",label=TRUE)
p3

p3 <- DimPlot(sample_m, reduction = "umap", split.by = "sex_dse_trt",label=TRUE)
p3
dev.off()

Idents(sample) <- "disease"
p3 <- DimPlot(sample, reduction = "umap", label=FALSE)
p3

Idents(sample) <- "trt"
p3 <- DimPlot(sample, reduction = "umap", label=FALSE)
p3

Idents(sample) <- "sex"
p3 <- DimPlot(sample, reduction = "umap", label=FALSE)
p3

Idents(sample) <- "disease_trt"
sample_ae <- subset(sample,idents=c("Asthma_E-cig"))
Idents(sample_ae) <- "sex"
levels(Idents(sample_ae))
# Reorder the Idents (for example, if you want to reorder them as 'B cells', 'T cells', 'NK cells')
new_order <- c("Female", "Male")
# Reassign the identities with the new order
Idents(sample_ae) <- factor(Idents(sample_ae), levels = new_order)
p3 <- DimPlot(sample_ae, reduction = "umap", label=FALSE)
p3

Idents(sample) <- "trt"
sample_ec <- subset(sample,idents=c("E-cig"))
Idents(sample_ec) <- "disease"
levels(Idents(sample_ec))
# Reorder the Idents (for example, if you want to reorder them as 'B cells', 'T cells', 'NK cells')
new_order <- c("Healthy", "Asthma")
# Reassign the identities with the new order
Idents(sample_ec) <- factor(Idents(sample_ec), levels = new_order)
p3 <- DimPlot(sample_ec, reduction = "umap", label=FALSE)
p3

Idents(sample) <- "disease"
sample_am <- subset(sample,idents=c("Asthma"))
Idents(sample_am) <- "trt"
levels(Idents(sample_am))
# Reorder the Idents (for example, if you want to reorder them as 'B cells', 'T cells', 'NK cells')
new_order <- c("Air", "E-cig")
# Reassign the identities with the new order
Idents(sample_am) <- factor(Idents(sample_am), levels = new_order)
p3 <- DimPlot(sample_am, reduction = "umap", label=FALSE)
p3


#find cluster
#marker genes 
DefaultAssay(sample) <- "RNA"
Idents(sample) <- "seurat_clusters"
clus_num <- unique(sample$seurat_clusters)
for (i in clus_num){
  
  marker <- paste0("marker", i)
  sample.marker <- FindMarkers(sample, ident.1=paste0(i),verbose = TRUE)
  
  write.csv(sample.marker,paste0("./cluster",marker,"_findmarkers.csv"), row.names = TRUE)
}