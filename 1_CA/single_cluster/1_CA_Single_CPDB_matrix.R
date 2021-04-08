# Symbol to ensembl
# setup
library("biomaRt");library("Seurat")
cpdb.outdir = '/media/yschang/W/cpdb/1_CA/'

# get data
Cancer.sbj <- readRDS(paste0(cpdb.outdir, 'Cancer.rds'))
# EPiCA.sbj <- readRDS('/home/yschang/CPDB/EPiCA.sbj')
# Cancer.sbj <- subset(EPiCA.sbj, idents = 'Cancer')
# load("/media/yschang/W/Merged/3_Split/Results/Split_Except_KBY_mito20.RData")


# 1. Ann.file from Seurat.obj
annotation.file <- paste0(cpdb.outdir, "Cancer_annotation.txt")
annotation <- data.frame(Cell=rownames(Cancer.sbj@meta.data), cell_type=Cancer.sbj@active.ident)
annotation$Cell <- gsub("-1", ".1", annotation$Cell)
write.table(annotation, file = annotation.file, quote = F, sep = "\t", row.names = F, col.names = T)

# 2. get count dataset
row.count <-data.frame(Cancer.sbj@assays$RNA@counts)
count_matrix <- apply(row.count, 2, function(x) (x/sum(x))*10000)

## change symbolic row.names into Ensembl.id
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
Symbols <- row.names(count_matrix)
df_gene <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol"), values = Symbols, mart= mart)

## merging count_matrix with df_gene (row.name vs hgnc_symbol)
### change count_matrix row.name to 1st column
count_matrix <-cbind(hgnc_symbol=rownames(row.count), count_matrix)
rownames(count_matrix) <- 1:nrow(count_matrix)
count_matrix <-as.data.frame(count_matrix)

### merge
df_count_matrix <-merge(df_gene, count_matrix, by= 'hgnc_symbol')

### check and remove duplication
df_count_matrix$hgnc_symbol[duplicated(df_count_matrix$hgnc_symbol)]  
df_count_matrix = df_count_matrix[-which(duplicated(df_count_matrix$hgnc_symbol)), ]

### remove Symbol and set row.names
df_count_matrix <-df_count_matrix[,-1]
df_count_matrix <- data.frame(df_count_matrix, row.names = T)

### export
count_matrix.file <- paste0(cpdb.outdir, "Cancer_count_matrix.txt")
write.table(df_count_matrix, file = count_matrix.file, quote = F, sep = "\t", col.names = T, row.names = T)

### save
saveRDS(Cancer.sbj, file = '/media/yschang/W/cpdb/1_CA/results/Cancer.rds')

rm(list=ls())
gc();q()
