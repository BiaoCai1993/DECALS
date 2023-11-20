path="/home/bc758/research/decals/"
library(Matrix)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

source(file=paste0(path,"code/supp/sc_analysis2.R"))
load(file=paste0(path,"data/sc/linear/data_list1.RData"))
load(file=paste0(path,"data/sc/linear/top_inds.RData"))
load(file=paste0(path,"data/sc/linear/ct_props.RData"))
load(file=paste0(path,"data/rosmap/sig.RData"))
load(file=paste0(path,"data/sc/linear/features.RData"))
Biao_genes <-rownames(sig)


major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
               'Astrocyte', 'Microglia','Endothelial')
plot_within_versus_between(Biao_genes, T, 'DECALS')

