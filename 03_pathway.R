

library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)


#load data
pbmc.data <- data.frame(fread("pbmc.data_standard.txt",header=T,sep="\t"),row.names=1)


#check species
msigdbr_show_species()
h <- msigdbr(species = "Homo sapiens",
             category = "C2")


h <- select(h, gs_name, gene_symbol) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) #或entrez_gene

#remove duplicate genes
gs <- lapply(h, unique)


#subset genesets
a<-names(gs)
gs<-gs[a]
gs <- gs[lapply(gs, length) > 0]
head(gs)


#GSVA
gsym.expr.1<-pbmc.data
gsva_es <- gsva(as.matrix(gsym.expr.1), gs, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=10)
gsva_es[1:3,1:3]
gsva_es_forshiny<-t(gsva_es)
#row.names(gsva_es_forshiny)<-paste("cell",row.names(gsva_es_forshiny),sep="")
write.table(gsva_es_forshiny,"pathway_all.txt",col.names=T,row.names=T,sep="\t",quote=F)
write.table((colnames(gsva_es_forshiny)),"pathway_all_names.txt",col.names=F,row.names=F,sep="\t",quote=F)


