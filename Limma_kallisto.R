setwd("/home/rudra/Research/NorthCott_MDB/Kallisto_output")

library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsembl(version = 102)
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'hsapiens_gene_ensembl',
                      version = 103)
attributes = listAttributes(ensembl)
attributes[1:5,]
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = ensembl)
#write.csv(t2g, "~/DataBase/Homo_sapiens/TranscriptVSgene_103.csv")
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
Tx <- dplyr::select(t2g, "target_id", "ext_gene")
Tx$ext_gene

Tx[Tx$ext_gene=="CASC15",]

list <- list.files("./")
list <- list[3:163]
path <- file.path(list, "abundance.tsv") 
all(file.exists(path)) 
sample <- read.csv("../metadata.csv")
sample <- as_tibble(sample)
sample$PAT.ID
coldata <- c()
for (i in list) {
  if (i %in% sample$PAT.ID) {
    print(i)
    coldata<- append(coldata,i)
  }
}

keep <- sample$PAT.ID %in% coldata
targets <- sample[keep,]
targets$PAT.ID
list
## import from kallisto
#t2g <- read.table("~/DataBase/Homo_sapiens/tx2gene.txt")
#t2g <- t2g[,c("V2", "V3")]
#colnames(t2g) <- c("target_id", "gene_name")


Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

### data wrangling
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
counts <- Txi_gene$counts
dim(counts)

colSums(counts)
colnames(counts) = list
write.csv(counts, "rawCount.csv")
#################################
#circ_df <- read.csv("~/Research/CIRCULAR_RNA/Northcott/Analysis/voom_count.csv", row.names = 1)
#col <- colnames(circ_df)

#drop <- setdiff(list, col)
#dim(counts)
#counts <- counts[,!(colnames(counts) %in% drop)]
#dim(counts)
##########################################

myDGEList <- DGEList(counts)
cpm <- edgeR::cpm(myDGEList)
dim(myDGEList)
keepers <- rowSums(cpm>5)>=1 #user defined
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
factor <- myDGEList.filtered.norm$samples
#write.csv(factor, "HostGeneTMMfactor.csv")
targets <- targets[targets$PAT.ID %in% colnames(counts),]

class <- as.factor(targets$Group)
design <- model.matrix(~0 + class)
colnames(design) <- levels(class)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
norm <- v.DEGList.filtered.norm$E
write.csv(norm, "voomNormalized.csv")
we <-grep("^LINC", rownames(norm))
norm[we,]
dim(norm)
top58<- read.csv('/home/rudra/Research/CIRCULAR_RNA/Northcott/Analysis/top62_with_geneName.csv')
top_genes <- top58$gene_name
top_genes <- unique(top_genes)[1:47]
#colnames(norm) = list

common <- c()
for (i in top_genes) {
  if (i %in% rownames(norm)) {
    print(i)
    common <- append(common,i)
  }
}

mt <- c("MTATP6P1","MTND2P28","MTRNR2L12")
norm <- counts[mt,]
norm <- t(norm)
norm
targets <- as.data.frame(targets)
rownames(targets) <- targets$PAT.ID
data <- cbind(targets, norm)
colnames(data)
data <- data[,c("Group"    , "MTATP6P1"  ,"MTND2P28" , "MTRNR2L12")]
## barplots
###


data$Group <- as.character(data$Group)
data$Group <- factor(data$Group, level =c("WNT", "SHH", "Group3", "Group4"))
setwd("/home/rudra/Desktop/")

for (i in mt ) {
  df <- data[,c("Group",i)]
  colnames(df) <- c("Group", "Gene")
  head(df)
  ggplot(df, aes(x=Group, y=Gene)) + theme_classic() +
    #geom_boxplot(fill = c("blue", "blue","green", "green","green", "green", "purple","purple", 
    #"purple", "red", "red", "red"))+
    geom_boxplot(fill = c("green","blue","orange","red")) +
    ggtitle(i)+
    theme(axis.text.x = element_text(size = 15, face = "bold", angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 15, face = "bold", angle = 0, vjust = 0, hjust=0)
          ,axis.title.x=element_text(size=15,face="bold", vjust = 0.5 )
          ,axis.title.y=element_text(size=15,face="bold", hjust = 0.5, vjust = 1.5 ),
          plot.title = element_text(size = 20, face = "bold"))+
    stat_summary(fun=mean, geom="point", shape=23, size=4) +
    xlab("Cancer Group") +
    ylab("Raw Count") +
    #theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(i,"Raw.pdf"), scale = 0.9)   
}


