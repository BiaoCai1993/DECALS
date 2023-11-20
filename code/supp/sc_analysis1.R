.libPaths("/home/bc758/R/x86_64-pc-linux-gnu-library/4.0")
library(Matrix)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/pseudo_bulk'
list.files(data_dir)

annot <- read.csv('/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/cell-annotation.csv')
ct_props <- table(annot$cell.type) / nrow(annot)
ct_props * 100

# Focus on four main cell types for now
major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
               'Astrocyte', 'Microglia','Endothelial')

subjects1 <- read.table(sprintf('%s/%s_subjects.txt', data_dir, major_cts[1]))[[1]]
subjects2 <- read.table(sprintf('%s/%s_subjects.txt', data_dir, major_cts[2]))[[1]]
all(subjects1 == subjects2)


load_data <- function(ct){
  # counts
  if(ct=="Excitatory Neurons"){
    counts1 <- readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
    counts2 <- readMM(sprintf('%s/%s_counts.txt', data_dir, "Inhibitory Neurons"))
    counts=counts1+counts2
  }else{
    counts <- readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
  }
  
  features <- read.table(sprintf('%s/%s_features.txt', data_dir, ct))[[1]]
  subjects <- read.table(sprintf('%s/%s_subjects.txt', data_dir, ct))[[1]]
  rownames(counts) <- features
  colnames(counts) <- subjects
  # sparisty
  sparsity <- Matrix::rowSums(counts == 0) / ncol(counts)
  # seq depth
  seq_depth <- Matrix::colSums(counts)
  # remove samples with no reads
  counts <- counts[, - which(seq_depth == 0)] # same for all four cell types
  seq_depth <- seq_depth[- which(seq_depth == 0)]
  # log normalized data
  dense_counts <- as.matrix(counts)
  scaled_data <- t(t(dense_counts) / seq_depth) * median(seq_depth)
  log_nor_data <- log(scaled_data + 1)
  # mean expression levels
  mean_exp <- rowMeans(log_nor_data)
  mean_exp_order <- length(mean_exp) - rank(mean_exp)
  return(list(counts = counts, log_nor_data = log_nor_data, scaled_data = scaled_data,
              sample_summary = list(seq_depth = seq_depth),
              gene_summary = list(sparsity = sparsity, mean_exp = mean_exp, mean_exp_order = mean_exp_order)))
}

data_list <- lapply(major_cts, function(ct) load_data(ct))
names(data_list) <- major_cts

dim(data_list[[1]]$counts)

sapply(data_list, function(dl) summary(dl$gene_summary$sparsity))
sapply(data_list, function(dl) summary(dl$sample_summary$seq_depth))


par(mfrow = c(2,2))
for(ct in major_cts){
  for(i in c(100, 1000, 10000, 20000)){
    x <- data_list[[ct]]$scaled_data[which(data_list[[ct]]$gene_summary$mean_exp_order == i), ]
    hist(x, main = sprintf('%s, %i', ct, i))
  }
}

output_prefix <- 'output/GO'
load_genes <- function(geneset){
  # load relevant data
  if(geneset != 'AD_risk_genes'){
    genes <- readRDS(sprintf('%s/%s_genes.rds', output_prefix, geneset))
  }else{
    genes <- readRDS(sprintf('%s/AD_genes_anlyzed.rds', 'output/AD_genes/multiplex_new'))
  }
  return(genes)
}
features <- read.table(sprintf('%s/%s_features.txt', data_dir, major_cts[1]))[[1]]

get_gene_set_rank <- function(genes, title = ''){
  rank_df <- lapply(major_cts, function(ct) data_list[[ct]]$gene_summary$mean_exp_order[match(genes, features)])
  rank_df_long <- do.call(c, rank_df)
  #colnames(rank_df) <- major_cts
  rank_df_long <- data.frame(rank = rank_df_long,
                             ct = rep(major_cts, each = length(genes)))
  g <- ggplot(rank_df_long) + geom_boxplot(aes(y = ct, x = rank, fill = ct)) +
    theme_classic(base_size= 18) +
    theme(legend.position = "none") +
    labs(title = title)
  print(g)
  return(rank_df_long)
}

AD_genes <- load_genes('AD_risk_genes')
Ex_genes <- load_genes('GOCC_EXCITATORY_SYNAPSE')
Oli_genes <- load_genes('GOCC_MYELIN_SHEATH')
Ast_genes <- load_genes('GOBP_ASTROCYTE_DIFFERENTIATION')
load('output/Biao/Biao_gene_set.RData')
Biao_genes <-rownames(sig)
print(head(Biao_genes))
gene_list <- list(AD_genes, Ex_genes, Oli_genes, Ast_genes, Biao_genes)
names(gene_list) <- c('AD', 'Ex', 'Oli', 'Ast', 'Biao')

for(i in 1:length(gene_list)){
  get_gene_set_rank(gene_list[[i]], names(gene_list)[i])
}

aaa=get_gene_set_rank(Biao_genes, "DECALS")

top_inds <- list()
for(top_n in c(10000, 15000, 20000)){
  top_k_inds <- sapply(major_cts, function(ct){
    data_list[[ct]]$gene_summary$mean_exp_order < top_n
  })
  top_inds[[as.character(top_n)]] <- apply(top_k_inds, 1, function(x) all(x))
}


#######################################################################
############### Pseudo Bulk Construction ##############################
#######################################################################
library(Matrix)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/pseudo_bulk'
list.files(data_dir)

annot <- read.csv('/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/cell-annotation.csv')

major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
               'Astrocyte', 'Microglia')

subjects1 <- read.table(sprintf('%s/%s_subjects.txt', data_dir, major_cts[1]))[[1]]
subjects2 <- read.table(sprintf('%s/%s_subjects.txt', data_dir, major_cts[2]))[[1]]
all(subjects1 == subjects2)


load_data <- function(ct){
  # counts
  if(ct=="Excitatory Neurons"){
    counts1 <- readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
    counts2 <- readMM(sprintf('%s/%s_counts.txt', data_dir, "Inhibitory Neurons"))
    counts=counts1+counts2
  }else{
    counts <- readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
  }
  
  features <- read.table(sprintf('%s/%s_features.txt', data_dir, ct))[[1]]
  subjects <- read.table(sprintf('%s/%s_subjects.txt', data_dir, ct))[[1]]
  rownames(counts) <- features
  colnames(counts) <- subjects
  # sparisty
  sparsity <- Matrix::rowSums(counts == 0) / ncol(counts)
  # seq depth
  seq_depth <- Matrix::colSums(counts)
  # remove samples with no reads
  id_0=which(seq_depth == 0)
  if(length(id_0)>0){
    counts <- counts[, - id_0] # same for all four cell types
    seq_depth <- seq_depth[- id_0]
  }
  # log normalized data
  dense_counts <- as.matrix(counts)
  scaled_data <- t(t(dense_counts) / seq_depth) * median(seq_depth)
  log_nor_data <- log(scaled_data + 1)
  # mean expression levels
  mean_exp <- rowMeans(log_nor_data)
  mean_exp_order <- length(mean_exp) - rank(mean_exp)
  return(list(counts = counts, log_nor_data = log_nor_data, scaled_data = scaled_data,
              sample_summary = list(seq_depth = seq_depth),
              gene_summary = list(sparsity = sparsity, mean_exp = mean_exp, mean_exp_order = mean_exp_order)))
}

data_list <- lapply(major_cts, function(ct) load_data(ct))
names(data_list) <- major_cts

dim(data_list[[1]]$counts)

load("~/research/decals/data/ROSMAP/sig.RData")
Biao_genes <-rownames(sig)


subjects=colnames(data_list[[1]]$counts)
ct_props_all=matrix(0,length(subjects),5)
for(i in 1:length(subjects)){
  ct_tab=table(annot$cell.type[annot$individualID==subjects[i]])
  ct_props_all[i,1]=ct_tab[names(ct_tab)=="Excitatory Neurons"]+ct_tab[names(ct_tab)=="Inhibitory Neurons"]
  ct_props_all[i,2]=ct_tab[names(ct_tab)=="Oligodendrocytes"]
  ct_props_all[i,3]=ct_tab[names(ct_tab)=="Astrocyte"]
  ct_props_all[i,4]=ct_tab[names(ct_tab)=="Microglia"]
  ct_props_all[i,5]=ct_tab[names(ct_tab)=="Endothelial"]
  ct_props_all[i,]=ct_props_all[i,]/sum(ct_props_all[i,])
}
colnames(ct_props_all)=major_cts

pseudo_bulk=matrix(0,159,436)
for(k in 1:5){
  pseudo_bulk=pseudo_bulk+data_list[[k]]$scaled_data[Biao_genes,]
}

load(file=paste0(path,"data/sc/linear/pseudo_bulk.RData"))
load(file=paste0(path,"data/sc/linear/ct_props_all1.RData"))
load(file=paste0(path,"data/sc/linear/features.RData"))
load(file=paste0(path,"data/sc/linear/sel_gene.RData"))
major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
               'Astrocyte', 'Microglia')

counts1 <- readMM(sprintf('%s/%s_counts.txt', data_dir, "Excitatory Neurons"))
counts2 <- readMM(sprintf('%s/%s_counts.txt', data_dir, "Inhibitory Neurons"))
counts=counts1+counts2
for(ct in major_cts[-1]){
  counts <- counts+readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
}
seq_depth <- Matrix::colSums(counts)
dense_counts <- as.matrix(counts)
scaled_data <- t(t(dense_counts) / seq_depth) * median(seq_depth)
pseudo_bulk1=scaled_data[match(Biao_genes,features),]
save(pseudo_bulk1,file="/home/bc758/research/decals/data/ROSMAP/sc/new/pseudo_bulk1.RData")

par(mfrow=c(2,3))
for(k in 1:5){
  result=cbind(ct_props_all[,k],frac0[,k])
  colnames(result)=c("TRUE","EST")
  boxplot(result,main=cell_names[k])
}
par(mfrow=c(1,1))