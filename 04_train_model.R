library(data.table)


#mouse
scrna.count<-data.frame(fread("~/data/integration3/molecules.integrate_mRNA.txt",header=T,sep="\t"),row.names=1)
pbmc.ident <- data.frame(fread("./../matrix/mouse/pbmc.ident_low.txt",header=T,sep="\t"))
#human
#scrna.count<-data.frame(fread("~/data/integration4/molecules.integrate_mRNA.txt",header=T,sep="\t"),row.names=1)
#pbmc.ident <- data.frame(fread("./../matrix/human/pbmc.ident_low.txt",header=T,sep="\t"))

scrna.count<-scrna.count[,pbmc.ident$samples]
#create metadata new new new new
metadata.integrate<-matrix(ncol=2,nrow=ncol(scrna.count))
metadata.integrate[,1]<-colnames(scrna.count)
metadata.integrate[,2]<-substr(colnames(scrna.count), 1, 2)
colnames(metadata.integrate)<-c("cellID","tech")
metadata.integrate<-data.frame(metadata.integrate)
row.names(metadata.integrate)<-metadata.integrate[,1]
metadata.integrate$tech <- as.factor(metadata.integrate$tech)
levels(metadata.integrate$tech)


scrna.count<-t(scrna.count)
scrna.1.stage<-data.frame(scrna.count)

#normalize matrix
scrna.1.stage<-log10(scrna.1.stage+1)



########
#using gini
gini.index  <-function(x){ 
  x <- sort(x)  
  G <- sum(x * 1L:length(x)) 
  G <- 2 * G/sum(x) - (length(x) + 1L)
  G/length(x)
}
#using marker gene sets
load("/www/data/TCA/matrix/mouse/pbmc.markers_reference.rda")
#load("/www/data/TCA/matrix/human/pbmc.markers_sct.rda")
pbmc.markers<-subset(pbmc.markers, avg_logFC > 0.25)
pbmc.markers$cluster<-as.numeric(as.character(pbmc.markers$cluster))
pbmc.markers$cluster<-pbmc.markers$cluster+1
unique(pbmc.markers$cluster)

geneSets<- vector('list', max(pbmc.markers$cluster))
for (i in 1:max(pbmc.markers$cluster)){
  pbmc.markers_sub<-subset(pbmc.markers, pbmc.markers$cluster == i)
  geneSets[[i]]<-row.names(pbmc.markers_sub)
  genes_tmp<-geneSets[[i]]
  genes_tmp<-intersect(genes_tmp, colnames(scrna.1.stage))
  geneSets[[i]]<-genes_tmp
}

gini.sum.names<-names(geneSets)
gini.sum<-character(0)
for (i in 1:length(geneSets)){
  if (i %% 10 == 0) print(i)
  genes_tmp<-geneSets[[i]]
  genes_tmp<-intersect(genes_tmp, colnames(scrna.1.stage))
  scrna.1.stage_sub<-scrna.1.stage[, genes_tmp]
  
  gini.tmp <- apply(scrna.1.stage_sub, 1, gini.index)
  gini.sum<-data.frame(cbind(gini.sum, gini.tmp))
  gini.sum[,i]<-as.numeric(as.character(gini.sum[,i]))
}
colnames(gini.sum)<-gini.sum.names
gini.sum<-data.frame(gini.sum)
gini.sum[1:3,1:3]

for(i in 1:ncol(gini.sum)){
  gini.sum[is.na(gini.sum[,i]), i] <- mean(gini.sum[,i], na.rm = TRUE)
}


library(dplyr)
save(gini.sum, file = "./ML/gini.sum_human.rda")
load(file = "./ML/gini.sum_human.rda")

gini.sum <- tibble::rownames_to_column(gini.sum, "samples")

scrna_sub_join<-dplyr::inner_join(pbmc.ident,gini.sum,by="samples")
#row.names(scrna_sub_join)<-scrna_sub_join$samples
scrna_sub_join<-scrna_sub_join[,-1];scrna_sub_join<-scrna_sub_join[,-1];colnames(scrna_sub_join)[1]<-"identity" #cluster
scrna_sub_join<-data.frame(scrna_sub_join)
#scrna_sub_join<-round(scrna_sub_join,4)
scrna_sub_join[1:3,1:3]
apply(is.na(scrna_sub_join), 2, which)



########
#using markers

#select features
marker.stage<-data.frame(fread("./../matrix/pbmc.markers.csv",header=T,sep=",",quote=""),row.names=1)
marker_sub<-subset(marker.stage,marker.stage$p_val_adj<0.001 & marker.stage$avg_logFC > 1)
markers.string<-unique(marker_sub$gene)
features<-intersect(markers.string,colnames(scrna.1.stage))

load("/www/data/TCA/matrix/human/pbmc.markers_sct.rda")
pbmc.markers<-subset(pbmc.markers, avg_logFC > 7.5 & p_val_adj<0.001)
markers.string<-pbmc.markers$gene
markers.string<-unique(markers.string)
features<-intersect(markers.string,colnames(scrna.1.stage))


#features<-c("Pou5f2","RP24-341O14.2","RP23-66E21.1","Phf2","Ssna1","Srek1","Ankar","Dnah8","RP23-440F7.2","RP24-300N16.3","Hnrnpm","Rps27a","Rdh11","Tchp","Ccp110","RP23-448A11.15","RP23-412L13.2","Gmcl1","Gm26799")
library(dplyr)
scrna_sub <- scrna.1.stage %>%  dplyr::select(one_of(features))
scrna_sub <- tibble::rownames_to_column(scrna_sub, "samples")

scrna_sub_join<-inner_join(pbmc.ident,scrna_sub,by="samples")
#scrna_sub_join<-scrna_sub_join[,-1];scrna_sub_join<-scrna_sub_join[,-1] #identity
scrna_sub_join<-scrna_sub_join[,-1];scrna_sub_join<-scrna_sub_join[,-1];colnames(scrna_sub_join)[1]<-"identity" #cluster
scrna_sub_join<-data.frame(scrna_sub_join)




# random split samples into training and validation
scrna_sub_join$identity <- as.factor(scrna_sub_join$identity)
index <- sample(2,nrow(scrna_sub_join),replace = TRUE,prob=c(0.7,0.3))
traindata <- scrna_sub_join[index==1,]
testdata <- scrna_sub_join[index==2,]

#### SVM
library(e1071)
cats_svm_model <- svm(identity~.,data=traindata,type = "C")
cats_svm_model
#save(cats_svm_model, file = "./ML/human_classification_marker_svm.RData")
#save(cats_svm_model, file = "./ML/mouse_classification_marker_svm.RData")
#save(cats_svm_model, file = "./ML/human_classification_gini.RData")
#save(cats_svm_model, file = "./ML/mouse_classification_gini.RData")


#### random forest
library(randomForest)
cats_svm_model <- randomForest(identity~.,data=traindata,ntree=100)
plot(cats_svm_model)
cats_svm_model
#save(cats_svm_model, file = "./ML/human_classification_marker_forest.RData")
save(cats_svm_model, file = "./ML/mouse_classification_marker_forest.RData")


# training
cats_svm_model_pred_1 <- predict(cats_svm_model,traindata[,-1])
cats_table_1 <- table(pred=cats_svm_model_pred_1,true=traindata[,1])
accuracy_1 <- sum(diag(cats_table_1))/sum(cats_table_1)
accuracy_1


# validation
cats_svm_model_pred_2 <- predict(cats_svm_model,testdata[,-1])
cats_table_2 <- table(pred=cats_svm_model_pred_2,true=testdata[,1])
#cats_table_2
accuracy_2 <- sum(diag(cats_table_2))/sum(cats_table_2)
accuracy_2
