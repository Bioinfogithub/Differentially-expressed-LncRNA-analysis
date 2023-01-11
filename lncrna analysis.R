#load libraries
library(BiocManager)
#BiocManager::install ("edgeR", force = TRUE)
library(edgeR)       
library(limma)
library(Biobase)
library(DESeq2)
#load data
setwd("E:/mtech project papers/main data")

#reading files
f1<- read.csv("GSM5049806_v1.txt", head = T, sep = "\t")
f2<- read.csv("GSM5049807_v2.txt", head = T, sep = "\t")
f3<- read.csv("GSM5049808_v3.txt", head = T, sep = "\t")
f4<- read.csv("GSM5049809_s1.txt",head = T, sep = "\t")
f5<- read.csv("GSM5049810_s2.txt",head = T, sep = "\t")
f6<- read.csv("GSM5049811_s3.txt",head = T, sep = "\t")
head(f1)
dim(f1)

deg_data = data.frame(nrow = 17488, ncol = 7)#making dataframe with 88659 row and 17 column
dim(deg_data)
head(deg_data)
#remove columns 
#deg_data<- subset(deg_data, select = -c(X))
deg_data = data.frame(f1$gene_id, f1$FPKM, f2$FPKM, f3$FPKM,f4$FPKM,f5$FPKM,f6$FPKM)
head(deg_data)
dim(deg_data)
summary(deg_data)

#Rename multiple columns by name
 
colnames(deg_data) <- c("gene_id","control1","control2","control3","SRSF10-overexpression1","SRSF10-overexpression2","SRSF10-overexpression3")        
head(deg_data)
View(deg_data)
#write.csv(deg_data,file = "Processed_FPKM_data.csv",quote = F)

#converting first column as a row in original dataset so that ncol(deg_data) == nrow(meta)
deg_data<- data.frame(deg_data[,-1], row.names = deg_data[,1])
#Removing rows having all zeros ???
deg_data<- deg_data[rowSums(deg_data[])>0,]
head(deg_data)
View(deg_data)

#step 1: preparing the count data ......
condition <- factor(c(rep("control",3),rep("SRSF10-overexpression",3)))# setting condition

#read in samples info....
meta <- data.frame(condition)

#making sure the row names in metadata matches to colnames in countdata
# are they in same order?
#all(colnames(countData) == row.names(metaData))
ncol(deg_data) == nrow(meta)
deg_data = ceiling(deg_data) #converting some values in assay that are not integers

#step2: Construct DESEQDataSet Object......
dds <- DESeqDataSetFromMatrix(
  countData = deg_data,
  design = ~ condition,
  colData = meta)

normCounts<-rlog(dds,blind = FALSE)#regularized log' transformation

#Now we're ready to run DESEQ function
dds <- DESeq(dds)

#estimateSizeFactors
#This calculates the relative library depth of each sample 

#estimateDispersions
#estimates the dispersion of counts for each gene 

#nbinomWaldTest
#calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs

res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table

#Summary of differential gene expression
summary(res) #summary of results

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

#plotcounts
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(1,3))# plotting multiple graph in single plot 

plotCounts(dds, gene="G9606_24417", intgroup="condition")
plotCounts(dds, gene="G9606_734", intgroup="condition")
plotCounts(dds, gene="G9606_30998", intgroup="condition")
plotCounts(dds, gene="G9606_32448", intgroup="condition")
plotCounts(dds, gene="G9606_14276", intgroup="condition")
plotCounts(dds, gene="G9606_13768", intgroup="condition")