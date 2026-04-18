library(Seurat)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
data<-readRDS('/users/PAS2783/aphmao2024/PinHsun/single cell/ecig.rds')


new.cluster.ids <- c("Club", "Basal", "Basal","Goblet","Mesenchymal","Mesothelial", "ATII","Goblet","Basal","Ciliated","Basal","Club","Ionocytes","Tuft","Deuterosomal","Ciliated")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
##Umap style 1
umap_1<-DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) + NoLegend()+
  theme(
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 15)
  )
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/umap1_updated.pdf", 
       plot = umap_1, 
       width = 13, 
       height = 10,  
       units = "in")

##Umap style 2
umap2<-DimPlot(data, reduction = "umap", label = FALSE, pt.size = 0.5)+
  theme(
    axis.title = element_text(size = 20, face = "bold"), 
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 25)
  )
ggsave("/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/umap2_updated.pdf", 
       plot = umap2, 
       width = 13, 
       height = 10,  
       units = "in")

data$cell_type <- Idents(data)


data$celltype.disese <- paste(data$cell_type, data$disease, sep = "_")
data$celltype.esig <- paste(data$cell_type, data$trt, sep = "_")


df_trt <- data.frame(cell.type = data$cell_type, smoking = data$trt )

df_trt %>% group_by(smoking, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(smoking)%>%mutate(proportion=count/sum(count))

df_disease <- data.frame(cell.type = data$cell_type, disease = data$disease )

df_disease %>% group_by(disease, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(disease)%>%mutate(proportion=count/sum(count))

prop_disease = df_disease %>% group_by(disease, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(disease)%>%mutate(proportion=count/sum(count))

prop_disease = as.data.frame(prop_disease)
pie(prop_disease$proportion[prop_disease$disease=="Asthma"],prop_disease$cell.type[prop_disease$disease=="Asthma"],cex=.5)
pie(prop_disease$proportion[prop_disease$disease=="Healthy"],prop_disease$cell.type[prop_disease$disease=="Healthy"],cex=.5)

prop_trt = df_trt %>% group_by(smoking, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(smoking)%>%mutate(proportion=count/sum(count))
prop_trt=as.data.frame(prop_trt)
pie(prop_trt$proportion[prop_trt$smoking=="Air"],prop_trt$cell.type[prop_trt$smoking=="Air"],cex=.5)
pie(prop_trt$proportion[prop_trt$smoking=="E-cig"],prop_trt$cell.type[prop_trt$smoking=="E-cig"],cex=.5)


Idents(data) <- data$celltype.disese

for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data, ident.1 = paste0(cell.type, '_Asthma'), ident.2 = paste0(cell.type, '_Healthy'))
  write.csv(de, file = paste0(dir.out, cell.type,"_disease.csv"))
  
}

Idents(data) <- data$celltype.esig


for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data, ident.1 = paste0(cell.type, '_Air'), ident.2 = paste0(cell.type, '_E-cig'))
  write.csv(de, file = paste0(dir.out, cell.type,"_esig.csv"))
  
}



## Asthma Air vs. Asthma E-Cig

data_asthma = subset(x = data, subset = disease == "Asthma")
##Cell proportion comparison
data_asthma$cell_type <- Idents(data_asthma)
data_asthma$celltype.disese <- paste(data_asthma$cell_type, data_asthma$disease, sep = "_")
data_asthma$celltype.esig <- paste(data_asthma$cell_type, data_asthma$trt, sep = "_")


df_asthma_esig <- data.frame(cell.type = data_asthma$cell_type, trt = data_asthma$trt )

prop_asthma_esig<-df_asthma_esig %>% group_by(trt, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(trt)%>%mutate(proportion=count/sum(count))

prop_asthma_esig = as.data.frame(prop_asthma_esig)

prop_asthma_esig_air=prop_asthma_esig|>filter(trt=="Air")|>mutate(proportion =round(proportion,digits = 3)*100)

# Get the positions
prop_asthma_esig_air2 <- prop_asthma_esig_air %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))




ggplot(prop_asthma_esig_air, aes(x = "" , y = count, fill = cell.type)) +
  geom_bar(stat = "identity", width = 1)+
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = prop_asthma_esig_air2,
                   aes(y = pos, label = paste0(proportion, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void()+ labs(title="Asthma Air")


prop_asthma_esig_esig=prop_asthma_esig|>filter(trt=="E-cig")|>mutate(proportion =round(proportion,digits = 3)*100)

# Get the positions
prop_asthma_esig_ecig2 <- prop_asthma_esig_esig %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))




ggplot(prop_asthma_esig_esig, aes(x = "" , y = count, fill = cell.type)) +
  geom_bar(stat = "identity", width = 1)+
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = prop_asthma_esig_ecig2,
                   aes(y = pos, label = paste0(proportion, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void()+ labs(title="Asthma E-cig")


##Gene differential expression (Updated to ident1 = E-cig, ident2 = Air
Idents(data_asthma) <- data_asthma$celltype.esig

dir.out = '/users/PAS2783/aphmao2024/PinHsun/single cell/DE_0417_2026/'
CELLtype <- unique(data_asthma$cell_type)
for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data_asthma, ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  write.csv(de, file = paste0(dir.out, cell.type,"_trt_cond_asthma_0417.csv"))
  
}


## Healthy Air vs. Healthy E-Cig

data_healthy = subset(x = data, subset = disease == "Healthy")
##Cell proportion comparison
data_healthy$cell_type <- Idents(data_healthy)
data_healthy$celltype.disese <- paste(data_healthy$cell_type, data_healthy$disease, sep = "_")
data_healthy$celltype.esig <- paste(data_healthy$cell_type, data_healthy$trt, sep = "_")


df_healthy_esig <- data.frame(cell.type = data_healthy$cell_type, trt = data_healthy$trt )

prop_healthy_esig<-df_healthy_esig %>% group_by(trt, cell.type) %>% summarise(count=n(),.groups = 'drop')%>%group_by(trt)%>%mutate(proportion=count/sum(count))

prop_healthy_esig = as.data.frame(prop_healthy_esig)

prop_healthy_esig_air=prop_healthy_esig|>filter(trt=="Air")|>mutate(proportion =round(proportion,digits = 3)*100)

# Get the positions
prop_healthy_esig_air2 <- prop_healthy_esig_air %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))




ggplot(prop_healthy_esig_air, aes(x = "" , y = count, fill = cell.type)) +
  geom_bar(stat = "identity", width = 1)+
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = prop_healthy_esig_air2,
                   aes(y = pos, label = paste0(proportion, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void()+ labs(title="Healthy Air")


prop_healthy_esig_esig=prop_healthy_esig|>filter(trt=="E-cig")|>mutate(proportion =round(proportion,digits = 3)*100)

# Get the positions
prop_healthy_esig_ecig2 <- prop_healthy_esig_esig %>% 
  mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))




ggplot(prop_healthy_esig_esig, aes(x = "" , y = count, fill = cell.type)) +
  geom_bar(stat = "identity", width = 1)+
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = prop_healthy_esig_ecig2,
                   aes(y = pos, label = paste0(proportion, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void()+ labs(title="Healthy E-cig")


##Gene differential expression
Idents(data_healthy) <- data_healthy$celltype.esig

CELLtype <- unique(data_healthy$cell_type)
dir.out  <- '/users/PAS2783/aphmao2024/PinHsun/single cell/DE_0417_2026/'
for (i in 1:length(CELLtype)){
  cell.type <- CELLtype[i]
  de <- FindMarkers(data_healthy, ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  write.csv(de, file = paste0(dir.out, cell.type,"_trt_cond_healthy_0417.csv"))
  
}


##Heatmap
data.markers <- FindAllMarkers(data, only.pos = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

heatmap = DoHeatmap(data, features = top5$gene) + NoLegend()
heatmap=heatmap + theme(axis.text=element_text(size=2))
ggsave('heatmap.png',plot = heatmap)

##Visualize the marker genes and the cluster 

markers.to.plot <-c('FOXI1','CFTR','V-ATPase','ASCL2','POU2F3','CALCA','ASCL1','KRT5','KRT14','TP63','KRT5','KRT4','KRT14','SCGB1A1','FUT4','CDC20B','FOXN4','DEUP1','TPPP3','FOXJ1','DNAH5','MUC5AC','FOXJ1','SPDEF','MUC5B','MUC5AC','FOXA2','FOXA3')
markers.to.plot<- unique(markers.to.plot)

DotPlot(data, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
DotPlot(data, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "trt") +
  RotatedAxis()

table(data$disease)
table(data$trt)

# Plot split by disease
dot_plot_disease <- DotPlot(data, 
                            features = markers.to.plot, 
                            cols = c("blue", "red"), 
                            dot.scale = 8, 
                            split.by = "disease") + RotatedAxis()

# Plot split by treatment
dot_plot_trt <- DotPlot(data, 
                        features = markers.to.plot, 
                        cols = c("blue", "red"), 
                        dot.scale = 8, 
                        split.by = "trt") + RotatedAxis()

head(dot_plot_disease$data)
head(dot_plot_trt$data)



dot_plot_disease+scale_size_continuous(
  name = "Percent Expressed", 
  range = c(1, 8),
  breaks = seq(0, 100, by = 5),   # This sets the ticks you asked for
  limits = c(0, 100)              # Ensure full range is shown
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dot_plot_trt+scale_size_continuous(
  name = "Percent Expressed", 
  range = c(1, 8),
  breaks = seq(0, 100, by = 5),   # This sets the ticks you asked for
  limits = c(0, 100)              # Ensure full range is shown
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






DotPlot(data, 
        features = markers.to.plot, 
        group.by = "cell_type", 
        dot.scale = 8) +
  RotatedAxis() +
  scale_size_continuous(
    name = "Percent Expressed", 
    range = c(1, 8),
    breaks = seq(0, 100, by = 5),   # This sets the ticks you asked for
    limits = c(0, 100)              # Ensure full range is shown
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(data, 
        features = markers.to.plot, 
        group.by = "cell_type", 
        dot.scale = 8) +
  RotatedAxis() +
  scale_size_continuous(
    name = "Percent Expressed", 
    range = c(1, 8),
    breaks = seq(0, 100, by = 5),   # This sets the ticks you asked for
    limits = c(0, 100)              # Ensure full range is shown
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##Split by different conditions
split_data <- SplitObject(data, split.by = "disease")  # or "trt"


# Plot Healthy
dot_healthy <- DotPlot(split_data$Healthy, features = markers.to.plot, group.by = "cell_type") +
  ggtitle("Healthy") +
  RotatedAxis()

dot_healthy
# Plot Asthma
dot_asthma <- DotPlot(split_data$Asthma, features = markers.to.plot, group.by = "cell_type") +
  ggtitle("Asthma") +
  RotatedAxis()

split_trt <- SplitObject(data, split.by = "trt")  

# Plot Healthy
dot_air <- DotPlot(split_trt$Air, features = markers.to.plot, group.by = "cell_type") +
  ggtitle("Air") +
  RotatedAxis()

dot_ecig <- DotPlot(split_trt$`E-cig`, features = markers.to.plot, group.by = "cell_type") +
  ggtitle("E-cig") +
  RotatedAxis()




##Volcano plot
## Asthma Air vs. Asthma E-Cig
# Output directory for volcano plots
Idents(data_asthma) <- data_asthma$celltype.esig

dir.out <- "/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/Volcano_updated/"

CELLtype <- unique(data_asthma$cell_type)

library(ggrepel)
for (cell.type in CELLtype) {
  
  de<-FindMarkers(data_asthma,ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  
  de$gene <- rownames(de)
  de$sig <- abs(de$avg_log2FC) > 1 & de$p_val_adj < 0.05
  
  top_ten_sig <- head(de[de$sig, ], 10)
  
  xlim_val <- max(abs(de$avg_log2FC))
  
  p<-ggplot(de, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = sig)) +
    geom_text_repel(data =  top_ten_sig, aes(label = gene), size = 6, box.padding = 0.5,point.padding = 0.3,max.overlaps = Inf) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))+theme_classic(base_size = 18)
  
  # Save the plot
  ggsave(filename = paste0(dir.out, cell.type, "_volcano_plot_cond_Asthma.pdf"),
         plot = p, width = 10, height = 10)
}



## Healthy Air vs. Healthy E-Cig
Idents(data_healthy) <- data_healthy$celltype.esig

CELLtype <- unique(data_healthy$cell_type)

dir.out <- "/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_Figures_0416_2026/Volcano_updated/"


for (cell.type in CELLtype) {
  
  de<-FindMarkers(data_healthy,ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  
  de$gene <- rownames(de)
  de$sig <- abs(de$avg_log2FC) > 1 & de$p_val_adj < 0.05
  
  top_ten_sig <- head(de[de$sig, ], 10)
  
  xlim_val <- max(abs(de$avg_log2FC))
  
  p<-ggplot(de, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = sig)) +
    geom_text_repel(data =  top_ten_sig, aes(label = gene), size = 6, box.padding = 0.5,point.padding = 0.3,max.overlaps = Inf) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))+theme_classic(base_size = 18)
  
  # Save the plot
  ggsave(filename = paste0(dir.out, cell.type, "_volcano_plot_cond_Healthy.pdf"),
         plot = p, width = 10, height = 10)
}


##Gene list for pathway analysis

## Asthma Air vs. Asthma E-Cig
# Output directory for volcano plots
Idents(data_asthma) <- data_asthma$celltype.esig

dir.out <- "/users/PAS2783/aphmao2024/PinHsun/single cell/Gene_list_for_pathway/"

CELLtype <- unique(data_asthma$cell_type)

##Significant counts by cell population
sig_counts_asthma <- c()


for (cell.type in CELLtype) {
  
  de<-FindMarkers(data_asthma,ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  
  de$gene <- rownames(de)
  de$sig <- de$p_val_adj < 0.05
  ##Check number of sig
  n_sig<-sum(de$sig)
  sig_counts_asthma[cell.type] <- n_sig
  
  # pick genes
  de_sig <- de[de$sig, ]
  if (n_sig <= 800) {
    top_genes <- de_sig$gene
  } else {
    ord<-order(abs(de_sig$avg_log2FC), decreasing = TRUE)
    top_genes <- rownames(de_sig)[head(ord,800)]
  }
  
  gene_lists[[cell.type]] <- top_genes
  out_path <- file.path(dir.out, paste0(cell.type, "_gene_list_asthma.csv"))
  write.table(top_genes,file = out_path,row.names = FALSE,col.names = FALSE,quote = FALSE,sep = ",")
  
}



## Healthy Air vs. Healthy E-Cig
# Output directory for volcano plots
Idents(data_healthy) <- data_healthy$celltype.esig

dir.out <- "/users/PAS2783/aphmao2024/PinHsun/single cell/Gene_list_for_pathway/"

CELLtype <- unique(data_healthy$cell_type)

##Significant counts by cell population
sig_counts_healthy <- c()
for (cell.type in CELLtype) {
  
  de<-FindMarkers(data_healthy,ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  
  de$gene <- rownames(de)
  de$sig <- de$p_val_adj < 0.05
  ##Check number of sig
  n_sig<-sum(de$sig)
  sig_counts_healthy[cell.type] <- n_sig
  
  # pick genes
  de_sig <- de[de$sig, ]
  if (n_sig <= 800) {
    top_genes <- de_sig$gene
  } else {
    ord<-order(abs(de_sig$avg_log2FC), decreasing = TRUE)
    top_genes <- rownames(de_sig)[head(ord, 800)]
  }
  
  out_path <- file.path(dir.out, paste0(cell.type, "_gene_list_healthy.csv"))
  write.table(top_genes,file = out_path,row.names = FALSE,col.names = FALSE,quote = FALSE,sep = ",")
  
}






Idents(data_healthy) <- data_healthy$celltype.esig

CELLtype <- unique(data_healthy$cell_type)

dir.out <- "/users/PAS2783/aphmao2024/PinHsun/single cell/Volcano_Plots/"


for (cell.type in CELLtype) {
  
  de<-FindMarkers(data_healthy,ident.1 = paste0(cell.type, '_E-cig'), ident.2 = paste0(cell.type, '_Air'))
  
  de$gene <- rownames(de)
  de$sig <- abs(de$avg_log2FC) > 1 & de$p_val_adj < 0.05
  
  top_ten_sig <- head(de[de$sig, ], 10)
  
  xlim_val <- max(abs(de$avg_log2FC))
  
  p<-ggplot(de, aes(x = avg_log2FC, y = -log10(p_val))) +
    geom_point(aes(color = sig)) +
    geom_text(data =  top_ten_sig, aes(label = gene), vjust = -0.5, size = 2)+
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))
  
  # Save the plot
  ggsave(filename = paste0(dir.out, cell.type, "_volcano_plot_cond_Healthy.png"),
         plot = p, width = 6, height = 5)
}


##Checking
top_list<-list()
ord_list<-list()
top_gene<-list()
asthma_dir<-'/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_DE_1021_2025/Asthma/'
for (i in CELLtype) {
  de_gene = read.csv(paste0(asthma_dir, i, "_trt_cond_asthma_updated.csv"))
  top_list[[i]]<-de_gene%>%arrange(desc(abs(avg_log2FC)))%>%head(10)
  ord_list[[i]] <- all(diff(de_gene$p_val_adj) >= 0)
  top_gene[[i]]<-de_gene%>%filter(abs(de_gene$avg_log2FC) > 1 & de_gene$p_val_adj < 0.05)%>%head(10)%>%select(X)
}
names(top_list)<- as.vector(CELLtype)

top_list<-list()
ord_list<-list()
top_gene<-list()
healthy_dir<-'/users/PAS2783/aphmao2024/PinHsun/single cell/Updated_DE_1021_2025/Healthy/'
for (i in CELLtype) {
  de_gene = read.csv(paste0(healthy_dir, i, "_trt_cond_healthy_updated.csv"))
  top_list[[i]]<-de_gene%>%arrange(desc(abs(avg_log2FC)))%>%head(10)
  ord_list[[i]] <- all(diff(de_gene$p_val_adj) >= 0)
  top_gene[[i]]<-de_gene%>%filter(abs(de_gene$avg_log2FC) > 1 & de_gene$p_val_adj < 0.05)%>%head(10)%>%select(X)
}