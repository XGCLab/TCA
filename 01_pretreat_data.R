#load packages
library(Seurat)
library(plyr)
library(dplyr)
library(data.table)
library(R.utils)

#read data
molecules.sperm_sd_GSE104556<-
  Read10X(data.dir = "./sperm_sd_GSE104556/out/")
colnames(molecules.sperm_sd_GSE104556)<-paste("d1-",colnames(molecules.sperm_sd_GSE104556),sep="")
row.names(molecules.sperm_sd_GSE104556)[1:3]
col1<-colnames(molecules.sperm_sd_GSE104556)
row1<-as.character(row.names(molecules.sperm_sd_GSE104556))

molecules.sperm_pg_GSE121904<-data.frame(fread("./sperm_pg_GSE121904/out/GSE121904_testis_maturation_scRNAseq_raw_mapped_counts.csv",header=T,sep=","),row.names=1)
colnames(molecules.sperm_pg_GSE121904)<-paste("d2-",colnames(molecules.sperm_pg_GSE121904),sep="")
row.names(molecules.sperm_pg_GSE121904)[1:3]
col2<-colnames(molecules.sperm_pg_GSE121904)
row2<-as.character(row.names(molecules.sperm_pg_GSE121904))

molecules.sperm_edc_GSE123119<-Read10X(data.dir = "./sperm_edc_GSE123119/GSE113293/outs/filtered_feature_bc_matrix")
colnames(molecules.sperm_edc_GSE123119)<-paste("d3-",colnames(molecules.sperm_edc_GSE123119),sep="")
row.names(molecules.sperm_edc_GSE123119)[1:3]
col3<-colnames(molecules.sperm_edc_GSE123119)
row3<-as.character(row.names(molecules.sperm_edc_GSE123119))

molecules.sperm_nc_GSE124904<-
  read.table("./sperm_nc_GSE124904/out/GSE124904_aggregate_gene_cell_matrix.txt"
             ,header=T,row.names=1)
colnames(molecules.sperm_nc_GSE124904)<-paste("d4-",colnames(molecules.sperm_nc_GSE124904),sep="")
row.names(molecules.sperm_nc_GSE124904)[1:3]
col4<-colnames(molecules.sperm_nc_GSE124904)
row4<-as.character(row.names(molecules.sperm_nc_GSE124904))

molecules.sperm_devcell_GSE112393<-
  data.table::fread("./sperm_devcell_GSE112393/out/GSE112393_MergedAdultMouseST25_DGE.txt"
             ,header=T,sep="\t")
molecules.sperm_devcell_GSE112393<-data.frame(molecules.sperm_devcell_GSE112393,row.names=1)
colnames(molecules.sperm_devcell_GSE112393)<-paste("d5-",colnames(molecules.sperm_devcell_GSE112393),sep="")
row.names(molecules.sperm_devcell_GSE112393)[1:3]
col5<-colnames(molecules.sperm_devcell_GSE112393)
row5<-as.character(row.names(molecules.sperm_devcell_GSE112393))

molecules.sperm_cr_GSE107644<-
  data.table::fread("/home/wuyingcheng/data/sperm_cr_GSE107644/rawcount/CR_UMIcount.txt"
                    ,header=T,sep="\t")
molecules.sperm_cr_GSE107644<-data.frame(molecules.sperm_cr_GSE107644,row.names=1)
colnames(molecules.sperm_cr_GSE107644)<-paste("d6-",colnames(molecules.sperm_cr_GSE107644),sep="")
row.names(molecules.sperm_cr_GSE107644)[1:3]
col6<-colnames(molecules.sperm_cr_GSE107644)
row6<-as.character(row.names(molecules.sperm_cr_GSE107644))


molecules.sperm_cell_GSE108097<-
  data.table::fread("./sperm_cell_GSE108097/out/Figure2-batch-removed-testis.txt"
                    ,header=T,sep="\t")
molecules.sperm_cell_GSE108097<-data.frame(molecules.sperm_cell_GSE108097,row.names=1)
colnames(molecules.sperm_cell_GSE108097)<-paste("d7-",colnames(molecules.sperm_cell_GSE108097),sep="")
row.names(molecules.sperm_cell_GSE108097)[1:3]
col7<-colnames(molecules.sperm_cell_GSE108097)
row7<-as.character(row.names(molecules.sperm_cell_GSE108097))

molecules.sperm_elife_GSE113293<-
  data.table::fread("./sperm_elife_GSE113293/out/raw/merge_count.txt"
                    ,header=T,sep="\t")
molecules.sperm_elife_GSE113293<-data.frame(molecules.sperm_elife_GSE113293,row.names=1)
colnames(molecules.sperm_elife_GSE113293)<-paste("d8-",colnames(molecules.sperm_elife_GSE113293),sep="")
row.names(molecules.sperm_elife_GSE113293)[1:3]
col8<-colnames(molecules.sperm_elife_GSE113293)
row8<-as.character(row.names(molecules.sperm_elife_GSE113293))

molecules.sperm_elife_GSE116001<-fread("./sperm_elife_GSE116001/GSE116001_MS17_WT_counts_merged.csv.gz",header=T,sep=",")
molecules.sperm_elife_GSE116001<-data.frame(molecules.sperm_elife_GSE116001,row.names=1)
colnames(molecules.sperm_elife_GSE116001)<-paste("d9-",colnames(molecules.sperm_elife_GSE116001),sep="")
row.names(molecules.sperm_elife_GSE116001)[1:3]
col9<-colnames(molecules.sperm_elife_GSE116001)
row9<-as.character(row.names(molecules.sperm_elife_GSE116001))




#combine into a large matrix
molecules.sperm_sd_GSE104556<-data.frame(molecules.sperm_sd_GSE104556)
molecules.sperm_pg_GSE121904<-data.frame(molecules.sperm_pg_GSE121904)
#molecules.sperm_edc_GSE123119<-data.frame(molecules.sperm_edc_GSE123119)
molecules.sperm_nc_GSE124904<-data.frame(molecules.sperm_nc_GSE124904)
molecules.sperm_devcell_GSE112393<-data.frame(molecules.sperm_devcell_GSE112393)
molecules.sperm_cr_GSE107644<-data.frame(molecules.sperm_cr_GSE107644)
molecules.sperm_cell_GSE108097<-data.frame(molecules.sperm_cell_GSE108097)
molecules.sperm_elife_GSE113293<-data.frame(molecules.sperm_elife_GSE113293)
molecules.sperm_elife_GSE116001<-data.frame(molecules.sperm_elife_GSE116001)

molecules.sperm_sd_GSE104556.name<-tibble::rownames_to_column(molecules.sperm_sd_GSE104556, "Genes")
molecules.sperm_pg_GSE121904.name<-tibble::rownames_to_column(molecules.sperm_pg_GSE121904, "Genes")
#molecules.sperm_edc_GSE123119.name<-tibble::rownames_to_column(molecules.sperm_edc_GSE123119, "Genes")
molecules.sperm_nc_GSE124904.name<-tibble::rownames_to_column(molecules.sperm_nc_GSE124904, "Genes")
molecules.sperm_devcell_GSE112393.name<-tibble::rownames_to_column(molecules.sperm_devcell_GSE112393, "Genes")
molecules.sperm_cr_GSE107644.name<-tibble::rownames_to_column(molecules.sperm_cr_GSE107644, "Genes")
molecules.sperm_cell_GSE108097.name<-tibble::rownames_to_column(molecules.sperm_cell_GSE108097, "Genes")
molecules.sperm_elife_GSE113293.name<-tibble::rownames_to_column(molecules.sperm_elife_GSE113293, "Genes")
molecules.sperm_elife_GSE116001.name<-tibble::rownames_to_column(molecules.sperm_elife_GSE116001, "Genes")



#filter mRNA
genes<-read.table("mart_export_mm9.txt",header=T,sep="\t")
genes_sub<-subset(genes,Gene.type %in% c("protein_coding", "Mt_tRNA"))


molecules.sperm_sd_GSE104556.name<-subset(molecules.sperm_sd_GSE104556.name,Genes %in% genes_sub$Gene.name)
molecules.sperm_pg_GSE121904.name<-subset(molecules.sperm_pg_GSE121904.name,Genes %in% genes_sub$Gene.name)
#molecules.sperm_edc_GSE123119.name<-subset(molecules.sperm_edc_GSE123119.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_nc_GSE124904.name<-subset(molecules.sperm_nc_GSE124904.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_devcell_GSE112393.name<-subset(molecules.sperm_devcell_GSE112393.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_cr_GSE107644.name<-subset(molecules.sperm_cr_GSE107644.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_cell_GSE108097.name<-subset(molecules.sperm_cell_GSE108097.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_elife_GSE113293.name<-subset(molecules.sperm_elife_GSE113293.name, Genes %in% genes_sub$Gene.name)
molecules.sperm_elife_GSE116001.name<-subset(molecules.sperm_elife_GSE116001.name, Genes %in% genes_sub$Gene.name)

dim(molecules.sperm_sd_GSE104556.name)
dim(molecules.sperm_pg_GSE121904.name)
#dim(molecules.sperm_edc_GSE123119.name)
dim(molecules.sperm_nc_GSE124904.name)
dim(molecules.sperm_devcell_GSE112393.name)
dim(molecules.sperm_cr_GSE107644.name)
dim(molecules.sperm_cell_GSE108097.name)
dim(molecules.sperm_elife_GSE113293.name)
dim(molecules.sperm_elife_GSE116001.name)





molecules.integrate <- molecules.sperm_sd_GSE104556.name
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_pg_GSE121904.name,by="Genes",all=TRUE)
#molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_edc_GSE123119.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_nc_GSE124904.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_devcell_GSE112393.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_cr_GSE107644.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_cell_GSE108097.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_elife_GSE113293.name,by="Genes",all=TRUE)
molecules.integrate <- merge(x=molecules.integrate,y=molecules.sperm_elife_GSE116001.name,by="Genes",all=TRUE)


molecules.integrate<-data.frame(molecules.integrate)


# exclude NA
molecules.integrate[is.na(molecules.integrate)] <- 0


#output
write.table(molecules.integrate,"molecules.integrate_mRNA.txt",col.names=T,row.names=F,sep="\t",quote=F)


