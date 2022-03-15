# Tutorial for multifactor designing 
##https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
setwd("/Users/rudramani/work/TIGR4/htseq_result/GENE/")
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(dplyr)
## import data ####
count <- read.csv("gene_counts.csv", row.names = 1)
count
dim(count)
meta <- read.csv("colData.csv", row.names = 1)
dim(meta)

mean(rownames(meta)==colnames(count))


meta$class <- as.factor(meta$class)
meta$copper <- as.factor(meta$copper)
meta$strain <- as.factor(meta$strain)
meta$zinc <- as.factor(meta$zinc)


### Filter low count genes
dim(count)
keep <- rowSums(count>0)>2  ## at least 2 sample with positive count

f <- count[keep,]
dim(f)

head(meta)
table(meta$zinc)
## Create Deseq2 object ####
dds <- DESeqDataSetFromMatrix(countData =f ,
                              colData = meta,
                              design = ~ copper)

levels(dds$copper)
# if reference is not at first try to relevel
#Normalization for sequencing depth for PCA and eploratory analysis ####
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

## log transform
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
## Or we can use below function 
log.norm <- normTransform(dds)
head(assay(log.norm),1)
head(log.norm.counts,1)
rs <- rowSums(counts(dds))
par(mfrow= c(1,2))
boxplot(log2(counts(dds)[rs > 0,] + 1)) # not normalized
boxplot(log.norm.counts[rs > 0,]) # normalized

# Stabilizing count variance ####
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)


### PCA plot #####
plotPCA(log.norm, intgroup="copper")
plotPCA(rld, intgroup="copper")
plotPCA(vsd, intgroup=c("strain", "copper", "zinc"))





#We can make this plot even nicer using custom code from the ggplot2 library:
library(ggplot2)
(data <- plotPCA(vsd, intgroup=c("strain","copper", "zinc"), returnData=TRUE))
data2 = data[grep("WT", rownames(data)),]

(percentVar <- 100*round(attr(data2, "percentVar"),2))

ggplot(data2, aes(PC1,PC2, col=strain, shape = zinc)) + geom_point(size = 3) +
  ylab(paste0("PC2: ",percentVar[2], " % variance"))+
  xlab(paste0("PC1: ",percentVar[1], " % variance"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  coord_fixed(ratio = 0.5)
## Cluster Denddograms

## Only selectiong "WT"
vsd2 = vsd[,grep("WT", rownames(data))]
par(mfrow = c(1,1))
plot(hclust(dist(t(assay(vsd2)))), 
     labels=colData(vsd2)$strain,
     color_labels = vsd2$copper)


## Plot the coloring
library("dendextend")
dend <- as.dendrogram(hclust(dist(t(assay(vsd2)))))
labels_colors(dend)
#labels_colors(dend) <- vsd2$copper
plot(dend)

## Differential Analysis ####

design(dds) <- ~ copper + zinc + strain + copper:zinc + copper:strain + copper:zinc:strain
design(dds)
dds <- DESeq(dds)
# Above interaction doesnot work [gives error full model matrix is less than full rank] if a occurance of a  particula combination is zero
# to avoid this 
ml <- model.matrix(~ 0 + copper + zinc + strain + copper:zinc + copper:strain + copper:zinc:strain, meta)
colnames(ml)
ml_df = as.data.frame(unname(ml)) # some of last combinatinos 44,45 etc are zeros
idx <- which(colSums(ml_df)!=0)
ml <- ml[,idx]
design(dds) <- ml

dds <- DESeq(dds)
levels(dds$strain)

resultsNames(dds)
res <-results(dds)
res
res1 <-results(dds, contrast = list("copper50", "copper0") )
res1 <- as.data.frame(res1)
res1 <- res1[res1$pvalue <= 0.05,]

res1 <- res1[order(-res1$log2FoldChange),]
head(res1)
write.csv(res1, "copper50vs0.csv")


assay(dds, normalize = TRUE)
write.csv(assay(vsd), "vsd_normalized.csv")


#### Designed experimental plan #####

design(dds) <- ~ 0 + class
design(dds)
dds <- DESeq(dds)

resultsNames(dds)
res <-results(dds)

#res1 <-results(dds, contrast = list("copper50", "copper0") )

wt_group <- c("classWT1", "classWT2", "classWT3", "classWT4", "classWT5",
           "classWT6", "classWT7", "classWT8" , "classWT9" )

de <- function(group1, group2){
  res1 <-results(dds, contrast = list(group2, group1) )
  res1 <- as.data.frame(res1)
  res1 <- res1[res1$pvalue <= 0.05,]
  res1 <- res1[order(-res1$log2FoldChange),]
  return(res1)
  
}



for (i in 1:(length(wt_group)-1)){
  for (j  in 2:length(wt_group) ) {
    if (i <j){
          name = paste0(wt_group[j],"_vs_", wt_group[i])
          print(name)
          #compare =   de(wt_group[i], wt_group[j])
          #write.csv(compare, paste0("Analysis/", name, ".csv"))
}}}


## Metal Specific analysis ##
#group <- c("classWT1", "classCopY1", "classCupA1",  "classWT10","classWT12")
#group <- c("classWT3", "classCopY2",  "classWT11","classWT13")
#group <- c("classWT5", "classCopY3",  "classCupA3")
#group <- c("classWT9", "classCopY4")
group <- c("classWT2", "classCupA2")
for (i in 2:length(group)){
      name = paste0(group[i],"_vs_", group[1])
      print(name)
      compare =   de(group[1], group[i])
      write.csv(compare, paste0("Analysis/", name, ".csv"))
    }

