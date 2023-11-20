path="/home/bc758/research/decals/"
load(file=paste0(path,"data/rosmap/sig.RData"))
Biao_genes <-rownames(sig)
library(Matrix)

major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
               'Astrocyte', 'Microglia','Endothelial')

load(file=paste0(path,"data/sc/linear/data_list1.RData"))
data_list1=data_list
for(k in 1:5){
  data_list1[[k]]=data_list[[k]][[3]]
}
data_list=data_list1
features=rownames(data_list[[1]])
#annot<-read.csv(file=paste0(path,"data/sc/linear/cell-annotation.csv"))
load(file=paste0(path,"data/sc/linear/annot.RData"))
cor_all=vector("list",5)
for(kk in 1:5){
  cell_prop=matrix(0,ncol(data_list[[kk]]),5)
  for(i in 1:ncol(data_list[[kk]])){
    id0=which(annot$individualID==colnames(data_list[[kk]])[i])
    cell_num=table(annot$cell.type[id0])
    neu=cell_num["Excitatory Neurons"]+cell_num["Inhibitory Neurons"]
    mic=cell_num["Microglia"]
    ast=cell_num["Astrocyte"]
    oli=cell_num["Oligodendrocytes"]
    end=cell_num["Endothelial"]
    cell_num1=c(neu,mic,ast,oli,end)
    cell_prop[i,]=cell_num1#/sum(cell_num1)
  }
  id1=match(Biao_genes,features)
  id2=is.na(cell_prop[,1])
  cell_cor=matrix(0,159,5)
  for(k in 1:5){
    for(j in 1:159){
      cell_cor[j,k]=cor(data_list[[kk]][id1[j],!id2],cell_prop[!id2,k])
    }
  }
  cor_all[[kk]]=cell_cor
}

id1=match(Biao_genes,features)
id0=match(colnames(data_list[[1]]),colnames(data_list[[5]]))
data_list[[5]]=data_list[[5]][,id0]
cor_all_mean=matrix(0,159,10)
i=0
for(k1 in 1:4){
  for(k2 in (k1+1):5){
    i=i+1
    for(j in 1:159){
      cor_all_mean[j,i]=cor(data_list[[k1]][id1[j],],data_list[[k2]][id1[j],])
    }
    
  }
}

par(mfrow=c(2,3))
for(k in 1:5){
  cell_cor=cor_all[[k]]
  colnames(cell_cor)=c("Neu","Mic","Ast","Oli","End")
  cell_cor1=cell_cor[,c(1,4,2,3,5)]
  boxplot(cell_cor1,ylim=c(-0.58,0.58),main=colnames(cell_cor1)[k],outline=FALSE)
  lines(x=c(0,5.5),y=c(0,0),lty=2,col=2)
}
#boxplot(cor_all_mean[1:1590],outline=FALSE,main="Cor(u_k1,u_k2)",ylim=c(-0.6,0.7))
#lines(x=c(0,10.5),y=c(0,0),lty=2,col=2)
par(mfrow=c(1,1))

load(file=paste0(path,"data/sc/linear/cell_prop.RData"))
cell_prop1=cell_prop[,c(1,4,2,3,5)]
library(tidyr) 
mat <- as.data.frame(round(cor(cell_prop1), 2))
mat$var1 <- c("Neu","Oli","Mic","Ast","End")
data <- gather(mat, key = "var2", value = "corr", -var1)
data$var2<-rep(c("Neu","Oli","Mic","Ast","End"),each=5)

library(RColorBrewer)
library(ggplot2)

my_color <- brewer.pal(5, "Spectral")[5:1]

ggplot(data, aes(var1, var2, fill = corr)) +
  geom_point(aes(size = abs(corr)), shape = 21, colour = "black") +
  scale_fill_gradientn(colours = my_color) +
  scale_size_area(max_size = 15, guide = FALSE)+
  geom_text(aes(label = corr), size = 3, colour = "black", alpha = 0.7)+
  theme(
    axis.title = element_blank()
  )






