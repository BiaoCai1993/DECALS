.libPaths("/home/bc758/R/x86_64-pc-linux-gnu-library/4.0")
library(Matrix)
library(Seurat)
library(dplyr)
sc_obj <- readRDS('/gpfs/gibbs/pi/zhao/eddie_shared_data/SUMMER22/ROSMAP/GeneExpression/snRNAseqPFC_BA10/seurat_obj.rds')
sc_counts <- GetAssayData(object = sc_obj, slot = "counts") 
sc_counts <- as.matrix(sc_counts) 
sc_counts[1:5, 1:5]
sc_obj@meta.data$celltype %>% table()
sc_obj@meta.data$projid %>% table()

cell.type=sc_obj@meta.data$celltype
id_opc=which(cell.type=="Opc")
id_per=which(cell.type=="Per")
id_use=(1:length(cell.type))[-c(id_opc,id_per)]

cell.types=cell.type[id_use]
cell.types=as.character(cell.types)
cell.types[cell.types=="Ex"]="Neu"
cell.types[cell.types=="In"]="Neu"
cell.types=as.factor(cell.types)
subject.name=sc_obj@meta.data$projid[id_use]

rna_metadata = read.table("/gpfs/gibbs/pi/zhao/eddie_shared_data/SUMMER22/ROSMAP/GeneExpression/snRNAseqPFC_BA10/Metadata/snRNAseqPFC_BA10_biospecimen_metadata.csv", header=TRUE,sep=",")
match_id=match(subject.name,rna_metadata[,1])
subject.names=rna_metadata$individualID[match_id]
save(cell.types,file="/home/bc758/research/decals/data/ROSMAP/bisque/cell.types.RData")
save(subject.names,file="/home/bc758/research/decals/data/ROSMAP/bisque/subject.names.RData")

sc_rna_counts=t(sc_counts)[id_use,]
save(sc_rna_counts,file="/home/bc758/research/decals/data/ROSMAP/bisque/sc_rna_counts.RData")

matrix_to_expressionSet <- function(mat, ...){
  if(!is.matrix(mat)) warning(deparse(substitute(mat)), " is not a matrix") else {
    featureData <- rownames(mat)
    featureData <- as(as.data.frame(featureData), "AnnotatedDataFrame")
    rownames(featureData) <- rownames(mat)
    phenoData <- colnames(mat)
    phenoData <- as(as.data.frame(phenoData), "AnnotatedDataFrame")
    rownames(phenoData) <- colnames(mat)
    ExpressionSet(assayData = mat, phenoData = phenoData, featureData = featureData) }
}

load("~/research/decals/data/ROSMAP/bulk.RData")
load("~/research/decals/data/ROSMAP/estFrac_RNA_5CT.Robj")
maker_genes=rownames(bulk)
colnames(bulk)=rownames(estFrac_RNA_5CT)

load(file="/home/bc758/research/decals/data/ROSMAP/bisque/cell.types.RData")
load(file="/home/bc758/research/decals/data/ROSMAP/bisque/subject.names.RData")
load(file="/home/bc758/research/decals/data/ROSMAP/bisque/sc_rna_counts.RData")
#sc_rna_counts=t(sc_rna_counts)
#save(sc_rna_counts,file="/home/bc758/research/decals/data/ROSMAP/bisque/sc_rna_counts.RData")
load(file="/home/bc758/research/decals/data/ROSMAP/bisque/id_na.RData")
id_use=(1:ncol(sc_rna_counts))[-id_na]
library(Biobase)
#library(SingleCellExperiment)
library(BisqueRNA)


#bulk.eset=matrix_to_expressionSet(bulk)
#sc.eset=matrix_to_expressionSet(sc_rna_counts)

bulk.eset <- Biobase::ExpressionSet(assayData = bulk)

sample.ids=colnames(sc_rna_counts[,id_use])
individual.labels=subject.names[id_use]
cell.type.labels=cell.types[id_use]
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=individual.labels,
                       cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=sc_rna_counts[,id_use],
                                  phenoData=sc.pdata)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset=bulk.eset, sc.eset=sc.eset, markers=maker_genes)

save(res,file="/home/bc758/research/decals/data/ROSMAP/bisque/res.RData")

