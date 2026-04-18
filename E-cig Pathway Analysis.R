library(Seurat)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)

data<-readRDS('/users/PAS2783/aphmao2024/PinHsun/single cell/ecig.rds')


new.cluster.ids <- c("Club", "Basal", "Basal","Goblet","Mesenchymal","Mesothelial", "ATII","Goblet","Basal","Ciliated","Basal","Club","Ionocytes","Tuft","Deuterosomal","Ciliated")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)


data$cell_type <- Idents(data)
data$celltype.disese <- paste(data$cell_type, data$disease, sep = "_")
data$celltype.esig <- paste(data$cell_type, data$trt, sep = "_")

##Asthma
data_asthma = subset(x = data, subset = disease == "Asthma")
data_asthma$cell_type <- Idents(data_asthma)
data_asthma$celltype.disese <- paste(data_asthma$cell_type, data_asthma$disease, sep = "_")
data_asthma$celltype.esig <- paste(data_asthma$cell_type, data_asthma$trt, sep = "_")
Idents(data_asthma) <- data_asthma$celltype.esig


##Healthy
data_healthy = subset(x = data, subset = disease == "Healthy")
data_healthy$cell_type <- Idents(data_healthy)
data_healthy$celltype.disese <- paste(data_healthy$cell_type, data_healthy$disease, sep = "_")
data_healthy$celltype.esig <- paste(data_healthy$cell_type, data_healthy$trt, sep = "_")
Idents(data_healthy) <- data_healthy$celltype.esig


#####Pathway Visualization

## DE analysis with ident 1 = E_cig and ident 2 = Air
CELLtype <- unique(data_healthy$cell_type)
dir.out  <- '/users/PAS2783/aphmao2024/PinHsun/single cell/DE_1020/DE_Healthy/'
for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data_healthy, ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  write.csv(de, file = paste0(dir.out, cell.type,"_trt_cond_healthy.csv"))
  
}

dir.out = '/users/PAS2783/aphmao2024/PinHsun/single cell/DE_1020/DE_Asthma/'
CELLtype <- unique(data_asthma$cell_type)
for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data_asthma, ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  write.csv(de, file = paste0(dir.out, cell.type,"_trt_cond_asthma.csv"))
  
}


##Reference gene list
senescence<-read.csv('/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/Senescence_Hallmarks_collective.csv')
senmayo <- senescence$SenMayo
senmayo <- senmayo[senmayo != ""]
senmayo <- unique(senmayo)
senmayo<-toupper(senmayo)


##Universal genes for entire data
Cell_int<-c("ATII", "Basal", "Ciliated", "Club", "Goblet")

data_asthma_sub <- subset(data_asthma, subset = cell_type %in% Cell_int)
Idents(data_asthma_sub)<- data_asthma_sub$celltype.esig

##Asthma universal genes for the data of 5 cell of interest 
asthma_assay <- DefaultAssay(data_asthma_sub)
genes_all_asthma <- Features(GetAssay(data_asthma_sub, assay = asthma_assay))
length(genes_all_asthma)
##17097 for universe genes

asthma_dir <- "/users/PAS2783/aphmao2024/PinHsun/single cell/DE_1020/DE_Asthma/"

p_asthma<-c()
Cell_int<-c("ATII", "Basal", "Ciliated", "Club", "Goblet")

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(asthma_dir, cell.type, "_trt_cond_asthma.csv"))
  de_gene$cell_pop = cell.type
  gene_list<-as.character(de_gene$X)
  
  n.ti <- length(intersect(senmayo, gene_list))
  n.t <- length(senmayo)
  n.i <- length(gene_list)
  n.all <- 17097
  p_val<-0.5*dhyper(n.ti, n.t, n.all-n.t, n.i )+phyper(n.ti, n.t, n.all-n.t, n.i, lower.tail=FALSE)
  
  print(cell.type)
  print(n.i)
  print(n.ti)
  print(intersect(senmayo, gene_list))
  print(p_val)
  p_asthma[cell.type] <- p_val
  
  
}


##Healthy

data_healthy_sub <- subset(data_healthy, subset = cell_type %in% Cell_int)
Idents(data_healthy_sub)<- data_healthy_sub$celltype.esig


##Healthy universal genes for the data of 5 cell of interest 
healthy_assay <- DefaultAssay(data_healthy_sub)
genes_all_healthy <- Features(GetAssay(data_healthy_sub, assay = healthy_assay))
length(genes_all_healthy)
##17097 for universe genes

##Healthy directory
healthy_dir <- "/users/PAS2783/aphmao2024/PinHsun/single cell/DE_1020/DE_Healthy/"

p_healthy<-c()
Cell_int<-c("ATII", "Basal", "Ciliated", "Club", "Goblet")

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(healthy_dir, cell.type, "_trt_cond_healthy.csv"))
  de_gene$cell_pop = cell.type
  gene_list<-as.character(de_gene$X)
  
  n.ti <- length(intersect(senmayo, gene_list))
  n.t <- length(senmayo)
  n.i <- length(gene_list)
  n.all <- 17097
  p_val<-0.5*dhyper(n.ti, n.t, n.all-n.t, n.i)+phyper(n.ti, n.t, n.all-n.t, n.i, lower.tail=FALSE)
  
  print(cell.type)
  print(n.i)
  print(n.ti)
  print(intersect(senmayo, gene_list))
  print(p_val)
  p_healthy[cell.type] <- p_val
  
  
}

names(p_asthma)=Cell_int
names(p_healthy)=Cell_int

library(ggplot2)
df <- data.frame(
  cell_type = c(Cell_int, Cell_int),
  condition = c(rep("Healthy", 5),
                rep("Asthma",  5)),
  pval = c(as.numeric(p_healthy),
           as.numeric(p_asthma))
)

df$neglog10p <- -log10(df$pval)

library(writexl)
write_xlsx(df, path = "/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/enrichment_table.xlsx")

enc<-ggplot(df, aes(cell_type, neglog10p, fill = condition)) +geom_bar(position="dodge", stat="identity")+
  geom_hline(yintercept = -log10(0.05),linetype=2)+geom_hline(yintercept = -log10(0.05/10),linetype=1)+
  labs(x = "Cell population", y = "-log10(p-value)", title = "Enrichment by cell type")+theme_classic(base_size = 25)
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/enrichment.pdf", 
       plot = enc, 
       width = 12,    
       height = 8,    
       units = "in")

##Focus on healthy population only
df_healthy<-df%>%filter(condition=='Healthy')
enc_healthy<-ggplot(df_healthy, aes(cell_type, neglog10p, fill = condition)) +geom_bar(position="dodge", stat="identity")+
  geom_hline(yintercept = -log10(0.05),linetype=2)+geom_hline(yintercept = -log10(0.05/5),linetype=1)+
  labs(x = "Cell population", y = "-log10(p-value)", title = "Enrichment by cell type")

ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/enrich_healthy.png", plot = enc_healthy)

##Heatmap
overlap_asthma<-list()

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(asthma_dir, cell.type, "_trt_cond_asthma.csv"))
  de_gene$cell_pop = cell.type
  gene_list<-as.character(de_gene$X)
  
  overlap_asthma[[cell.type]] <-intersect(senmayo, gene_list)
  
  
  
}

genes_int_asthma <- Reduce(intersect, overlap_asthma)
asthma_int_genes <- sort(as.character(genes_int_asthma))

mat_padj_asthma <- matrix(0,
                          nrow = length(asthma_int_genes),
                          ncol = length(Cell_int),
                          dimnames = list(asthma_int_genes, Cell_int)
)

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(asthma_dir, cell.type, "_trt_cond_asthma.csv"))
  selected<-de_gene%>%filter(X%in%genes_int_asthma)%>%select(X,p_val_adj)
  selected<-selected%>%arrange(X)
  
  mat_padj_asthma[, cell.type] <- -log10(selected$p_val_adj)
}
mat_padj_asthma[ mat_padj_asthma > 30] <- 30
library(pheatmap)
pheatmap(mat_padj_asthma,
         main = "Asthma -log10(adjusted p-value)",
         cluster_rows = FALSE, cluster_cols = FALSE,
         na_col = "grey90")



##Healthy
overlap_healthy<-list()

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(healthy_dir, cell.type, "_trt_cond_healthy.csv"))
  de_gene$cell_pop = cell.type
  gene_list<-as.character(de_gene$X)
  
  overlap_healthy[[cell.type]] <-intersect(senmayo, gene_list)
  
  
  
}

genes_int_healthy <- Reduce(intersect, overlap_healthy)


genes <- sort(as.character(genes_int_healthy))
cts   <- Cell_int                         

mat_padj_healthy <- matrix(
  0,
  nrow = length(genes),
  ncol = length(cts),
  dimnames = list(genes, cts)
)

for (cell.type in Cell_int) {
  
  de_gene = read.csv(paste0(healthy_dir, cell.type, "_trt_cond_healthy.csv"))
  selected<-de_gene%>%filter(X%in%genes_int_healthy)%>%select(X,p_val_adj)
  selected<-selected%>%arrange(X)
  
  mat_padj_healthy[, cell.type] <- -log10(selected$p_val_adj)
}

mat_padj_healthy[ mat_padj_healthy > 30] <- 30
pheatmap(mat_padj_healthy,
         main = "Healthy -log10(adjusted p-value)",
         cluster_rows = FALSE, cluster_cols = FALSE,
         na_col = "grey90")



##Dots plot
##Subset data to focus on the 5 cell pop

data_sub<-subset(data, subset = cell_type %in% Cell_int)
Idents(data_sub)<- data_sub$celltype.esig
gene_itrs<-c('TFAM', 'CDKN2A', 'CDKN1A','CCL2', 'SERPINE1', 'IL6')
dp<-DotPlot(data_sub, features = gene_itrs, cols = c("blue", "red"), dot.scale = 8, split.by = "disease") +
  RotatedAxis()
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/dot_plot.png", plot = dp,
       width = 12, height = 10, units = "in", dpi = 300)


##Feature plot

fp1<-FeaturePlot(data_sub, features = gene_itrs[1:3], split.by = "disease", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/featureplots_1.png", plot = fp1,
       width = 12, height = 10, units = "in", dpi = 300)

fp2<-FeaturePlot(data_sub, features = gene_itrs[4:6], split.by = "disease", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/featureplots_2.png", plot = fp2,
       width = 12, height = 10, units = "in", dpi = 300)


##Adjust order of X axia in the violin plot 

Idents(data_sub) <- factor(Idents(data_sub), levels = c('ATII_Air','ATII_E-cig','Basal_Air','Basal_E-cig','Ciliated_Air','Ciliated_E-cig','Club_Air','Club_E-cig','Goblet_Air','Goblet_E-cig'))
vio_plots_1 <- VlnPlot(data_sub, features = gene_itrs[1:3], split.by = "disease",
                       pt.size = 0, combine = FALSE)
font_size <- 20
v1<-wrap_plots(plots = vio_plots_1, ncol = 1)

ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/violin_1.png", plot = v1,
       width = 12, height = 10, units = "in", dpi = 300)

##Updated on 4/26/2026
Idents(data_sub) <- factor(Idents(data_sub), levels = c('ATII_Air','ATII_E-cig','Basal_Air','Basal_E-cig','Ciliated_Air','Ciliated_E-cig','Club_Air','Club_E-cig','Goblet_Air','Goblet_E-cig'))
vio_plots_1 <- VlnPlot(data_sub, features = gene_itrs[2:3], split.by = "disease",
                       pt.size = 0, combine = FALSE)

font_size <- 20

v1 <- wrap_plots(plots = vio_plots_1, ncol = 1) & 
  theme_classic() & 
  theme(
    text = element_text(size = font_size),
    axis.title.y = element_text(size = font_size),
    axis.text.y = element_text(size = font_size - 4),
    axis.text.x = element_text(size = font_size - 2, angle = 45, hjust = 1),
    plot.title = element_text(size = font_size + 2, face = "bold"),
    legend.text = element_text(size = font_size - 4)
  )

ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/violin_updated.pdf", 
       plot = v1, 
       width = 10, 
       height = 15,  
       units = "in")








vio_plots_2 <- VlnPlot(data_sub, features = gene_itrs[4:6], split.by = "disease",
                       pt.size = 0, combine = FALSE)
v2<-wrap_plots(plots = vio_plots_2, ncol = 1)
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/violin_2.png", plot = v2,
       width = 12, height = 10, units = "in", dpi = 300)


##Violin plot focusing on healthy population
data_sub_healthy<-subset(data_sub, subset = disease == "Healthy")
Idents(data_sub_healthy) <- factor(Idents(data_sub_healthy), levels = c('ATII_Air','ATII_E-cig','Basal_Air','Basal_E-cig','Ciliated_Air','Ciliated_E-cig','Club_Air','Club_E-cig','Goblet_Air','Goblet_E-cig'))

vio_plots_1_healthy <- VlnPlot(data_sub_healthy,features = gene_itrs[2:3],pt.size = 0,combine = TRUE)+theme(legend.position = "none")
vio_plots_1_healthy  <- (vio_plots_1_healthy  + plot_layout(ncol = 1)) & theme(legend.position = "none")

ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/violin_1_healthy.pdf", plot = vio_plots_1_healthy,
       width = 12, height = 10, units = "in", dpi = 300)

vio_plots_2_healthy <- VlnPlot(data_sub_healthy,features = gene_itrs[4:6],pt.size = 0,combine = TRUE)+theme(legend.position = "none")
vio_plots_2_healthy  <- (vio_plots_2_healthy  + plot_layout(ncol = 1)) & theme(legend.position = "none")

ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/pathway_visualization/violin_2_healthy.png", plot = vio_plots_2_healthy,
       width = 12, height = 10, units = "in", dpi = 300)

