# set-up
rm(list = ls()); gc()
library(Seurat)
library(dplyr)
Cancer.outdir = '/media/yschang/W/cpdb/1_CA/results/'

# get previous work
load('/media/yschang/W/Merged/3_Split/Results/Split_Except_KBY_mito20.RData')
Cancer.sbj <- subset(YS_Itg.sfj, idents = 'NE')

# Define useful fxs
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

# color and markers
Cancer.color <- c('#CCCCCC', '#C0C0C0', '#F8766D')
Cancer.mk <- c("SPINK1", "NAPSA", "HOPX", "LPCAT1", "AREG", "LAPTM4B")

# Check Cancer identity
DimPlot(object = YS_Itg.sfj, cells.highlight = WhichCells(Cancer.sbj), cols.highlight = '#F8766D', reduction = "umap", label = T)
FeaturePlot(YS_Itg.sfj, cols = Cancer.color, features = Cancer.mk)

#################
# 1. Reclustering
#################
DefaultAssay(Cancer.sbj) <- 'RNA'
Cancer.sbj <- FindVariableFeatures(Cancer.sbj, selection.method = "vst", nfeatures = 2000, assay = 'RNA'); v.genes <- VariableFeatures(Cancer.sbj)

DefaultAssay(Cancer.sbj) <- 'integrated'
Cancer.ssj <- ScaleData(object = Cancer.sbj, features = rownames(Cancer.sbj)); rm(Cancer.sbj)
Cancer.spj <- RunPCA(object = Cancer.ssj, features = v.genes); rm(Cancer.ssj)

# Estimate dimension
ElbowPlot(object = Cancer.spj)
Cancer.spj <- JackStraw(Cancer.spj, num.replicate = 100); Cancer.spj <- ScoreJackStraw(Cancer.spj, dims = 1:20); JackStrawPlot(Cancer.spj, dims = 1:20)

# UMAP
Cancer.suj <- RunUMAP(Cancer.spj, dims = 1:5)# rm(Cancer.spj)

# Clustering
Cancer.suj <- FindNeighbors(Cancer.suj, dims = 1:10)
Cancer.suj <- FindClusters(Cancer.suj, resolution = 0.1)
DimPlot(Cancer.suj, reduction = "umap", label = T, group.by = 'orig.ident')
DimPlot(Cancer.suj, reduction = "umap", label = F, pt.size = 1.0)

# Find markers
DefaultAssay(Cancer.suj) <- 'RNA'

# FindAllmarkers
Cancer_all.mk <- FindAllMarkers(Cancer.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(Cancer_all.mk, file = paste0(Cancer.outdir, "Integ_Cancer_all.mk"), quote = F, row.names = FALSE, sep = "\t")

Cancer_all.pw <- FindAllMarkers(Cancer.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
Cancer_all <- merge(Cancer_all.mk, Cancer_all.pw, by="row.names", all.x = F, all.y = F)
Cancer_all <- subset(Cancer_all, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
Cancer_all <- Cancer_all[rev(order(Cancer_all$cluster.x, Cancer_all$avg_log2FC.x, Cancer_all$power)), ]
write.table(x = Cancer_all, file = paste0(Cancer.outdir, "/Integ_Cancer_all.mp"), quote = F, sep = "\t")

# visualize by heatmap
Cancer_Top3.mk <- Cancer_all.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(Cancer.suj, features = Cancer_Top3.mk$gene, slot = 'count') + NoLegend()

#####################
# 2. Assign active id
#####################

# Keep subcluster.id
Cancer.suj$subcluster.id <- Idents(Cancer.suj)

# Cancer.sfj
Cancer.sfj <- Cancer.suj
Cancer.sfj <- RenameIdents(Cancer.sfj, "0" = "CA_0", "1" = "CA_1",
                           "2" = "CA_2", "3" = "CA_3", "4" = "CA_4",
                           "5" = "CA_5", "6" = "CA_6")


PlotClusterTree(BuildClusterTree(Cancer.sfj))
DimPlot(Cancer.sfj, reduction = "umap")
rm(Cancer.suj); rm(Cancer.spj)
##############
# 3. Visualize
##############
# Original DimPlot
jpeg(filename = paste0(Cancer.outdir, "Cancer_Dimplot.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Cancer.sfj, reduction = "umap", label = F, pt.size = 0.1)
dev.off()

# Save as EPS
Plt <- DimPlot(Cancer.sfj, reduction = "umap", label = F, pt.size = 0.1)
postscript(paste0(Cancer.outdir, "Cancer_Dimplot.eps"))
Plt
dev.off()

################################################################################
### Save EndProduct ############################################################
################################################################################

save.image(paste0(Cancer.outdir, "1_Integ_Cancer.RData"))
saveRDS(Cancer.sfj, paste0(Cancer.outdir, "Cancer.rds"))

rm(list=ls())
gc();q()
