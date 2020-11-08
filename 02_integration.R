
library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)


#read data
molecules.integrate<-data.frame(fread("molecules.integrate_mRNA.txt",header=T,sep="\t"),row.names=1)

options(future.globals.maxSize = 40000 * 1024^2)



#create metadata
metadata.integrate<-matrix(ncol=2,nrow=ncol(molecules.integrate))
metadata.integrate[,1]<-colnames(molecules.integrate)
metadata.integrate[,2]<-substr(colnames(molecules.integrate), 1, 2)
colnames(metadata.integrate)<-c("cellID","tech")
metadata.integrate<-data.frame(metadata.integrate)
row.names(metadata.integrate)<-metadata.integrate[,1]
metadata.integrate$tech <- as.factor(metadata.integrate$tech)
levels(metadata.integrate$tech)

metadata.integrate<-subset(metadata.integrate, metadata.integrate$tech %in% c("d1", "d2", "d5", "d6", "d7", "d8", "d9"))
samples<-row.names(metadata.integrate)
molecules.integrate<-molecules.integrate[,samples]


#create Seurat object and integrate
obj <- CreateSeuratObject(molecules.integrate, min.cells = 3, min.features = 200, meta.data = metadata.integrate)
obj.list <- SplitObject(obj, split.by = "tech")

all_features <- lapply(obj.list, row.names) %>% Reduce(intersect, .)
save(all_features, file = "pbmc.all_features_reference.rda")


for (i in 1:length(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)
}

obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 5000)
print("SelectIntegrationFeatures")
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features,  verbose = FALSE)
print("PrepSCTIntegration")

reference_dataset <- which(names(obj.list) == "d5")


obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = obj.features, verbose = FALSE, reference = reference_dataset)
print("FindIntegrationAnchors")
obj.integrated <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", features.to.integrate = all_features)
print("IntegrateData")
dim(obj.integrated)
length(all_features)



DefaultAssay(obj.integrated) <- "integrated"


#QC
obj.integrated_tmp<-obj.integrated@assays$integrated@data
obj.integrated_tmp <- CreateSeuratObject(obj.integrated_tmp, min.cells = 0, min.features = 0)
obj.integrated_tmp[["percent.mt"]] <- PercentageFeatureSet(obj.integrated_tmp, pattern = "^Mt.")
obj.integrated_tmp[["percent.rps"]] <- PercentageFeatureSet(obj.integrated_tmp, pattern = "^Rps")
obj.integrated_tmp[["percent.rpl"]] <- PercentageFeatureSet(obj.integrated_tmp, pattern = "^Rpl")

obj.integrated[["percent.mt"]] <- obj.integrated_tmp[["percent.mt"]]
obj.integrated[["percent.rps"]] <- obj.integrated_tmp[["percent.rps"]]
obj.integrated[["percent.rpl"]] <- obj.integrated_tmp[["percent.rpl"]]



#obj.integrated[["percent.mt"]] <- PercentageFeatureSet(obj.integrated, pattern = "^Mt.")
#obj.integrated[["percent.rps"]] <- PercentageFeatureSet(obj.integrated, pattern = "^Rps")
#obj.integrated[["percent.rpl"]] <- PercentageFeatureSet(obj.integrated, pattern = "^Rpl")
#VlnPlot(obj.integrated, features = c("percent.mt", "percent.rps", "percent.rpl"), ncol = 3, group.by = "identity.x")

dim(obj.integrated)
#obj.integrated <- subset(obj.integrated, subset = percent.rps < 4 & percent.rpl < 4)
dim(obj.integrated)

all.genes <- rownames(obj.integrated)
#obj.integrated <- ScaleData(obj.integrated, features = all.genes, vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl"))
#obj.integrated <- ScaleData(obj.integrated, features = all.genes, vars.to.regress = c("percent.rps", "percent.rpl"))
obj.integrated <- ScaleData(obj.integrated, features = all.genes, vars.to.regress = c("percent.mt"))
#obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
print("scale")

obj.integrated <- FindVariableFeatures(object = obj.integrated)
obj.integrated <- RunPCA(obj.integrated, features = VariableFeatures(object = obj.integrated))
print("pca")

ElbowPlot(obj.integrated)

pca_to_use <- 20 
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:pca_to_use)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.5)
head(Idents(obj.integrated), 5)
print("FindClusters")


#dimention reduction
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:30)
print("RunUMAP")

obj.integrated <- RunTSNE(object = obj.integrated, dims = 1:30, check_duplicates = FALSE)
print("RunTSNE")


pbmc.markers <- FindAllMarkers(obj.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(pbmc.markers, file = "pbmc.markers_reference.rda")
print("FindAllMarkers")


pbmc.ident <- Idents(obj.integrated)
save(pbmc.ident, file = "pbmc.ident_reference.rda")
print("Idents")



#plot
pdf("obj.integrated.standard_reference.pdf",height=5,width=15)
#plot1 <- DimPlot(obj.integrated, reduction = "umap", group.by = "tech")
#plot2 <- DimPlot(obj.integrated, reduction = "pca", group.by = "tech")
#plot3 <- DimPlot(object = obj.integrated, reduction = "tsne",group.by = "tech")
#CombinePlots(plots = list(plot1, plot2, plot3))
#plot1
DimPlot(obj.integrated, reduction = "umap", label = TRUE, pt.size = .01,split.by = "tech")
DimPlot(obj.integrated, reduction = "umap", label = TRUE, pt.size = .01,group.by = "tech")
DimPlot(object = obj.integrated, reduction = "tsne",pt.size = .01,group.by = "tech")
DimPlot(object = obj.integrated, reduction = "tsne",pt.size = .01,split.by = "tech")
DimPlot(obj.integrated, reduction = "pca",pt.size = .01, group.by = "tech")

dev.off()
dev.off()





#output tsne, umao, and pca
tsne<-Embeddings(object = obj.integrated[["tsne"]])
umap<-Embeddings(object = obj.integrated[["umap"]])
pca<-Embeddings(object = obj.integrated[["pca"]])

write.table(tsne,"tsne.txt",col.names=T,row.names=T,sep="\t",quote=F)
write.table(pca[,1:2],"pca.txt",col.names=T,row.names=T,sep="\t",quote=F)
write.table(umap,"umap.txt",col.names=T,row.names=T,sep="\t",quote=F)




# Run the standard workflow for visualization and clustering
obj.integrated <- RunALRA(object = obj.integrated)
print("RunALRA")
dim(obj.integrated)

#output data
pbmc.exp<-obj.integrated@assays$alra@data;pbmc.exp<-data.frame(pbmc.exp)
pbmc.exp<-round(pbmc.exp,4)
pbmc.exp<-tibble::rownames_to_column(pbmc.exp, "Genes")

write.table(pbmc.exp,"pbmc.data_reference_alra.txt",col.names=T,row.names=F,sep="\t",quote=F)
print("output")

#output data
pbmc.exp<-obj.integrated@assays$integrated@data;pbmc.exp<-data.frame(pbmc.exp)
pbmc.exp<-round(pbmc.exp,4)
pbmc.exp<-tibble::rownames_to_column(pbmc.exp, "Genes")

write.table(pbmc.exp,"pbmc.data_reference_integrateddata.txt",col.names=T,row.names=F,sep="\t",quote=F)
print("output")


#output data
pbmc.exp<-obj.integrated@assays$integrated@scale.data;pbmc.exp<-data.frame(pbmc.exp)
pbmc.exp<-round(pbmc.exp,4)
pbmc.exp<-tibble::rownames_to_column(pbmc.exp, "Genes")

write.table(pbmc.exp,"pbmc.data_reference_integratedscaledata.txt",col.names=T,row.names=F,sep="\t",quote=F)
print("output")




