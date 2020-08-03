"Here we are choosing a breast cancer dataset downloaded 
from Firehose (https://gdac.broadinstitute.org/).
These are next generation sequncing data that are 
provided already normalized."


setwd("~/work/Bioinformatics/")

"This dataset has genes as rows and patients as columns.
Patient IDs refer to TCGA IDs such as TCGA.3C.AAAU.01A.11R.A41B.07, etc. 
The top two lines of the file are header information"

mrnaNorm <- read.table("DataSets/Breast_cancer/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                       header = F, fill = T, skip = 0)
head(mrnaNorm)
dim(mrnaNorm)
col(mrnaNorm)
mrnaNorm[1:5, 1:5]
#We also want the patient IDs in the first row, descard firt two columns
mrnaIDs <- mrnaNorm[1,-c(1,2)]

# remove the header information
mrnaNorm <- mrnaNorm[-c(1,2),]
mrnaNorm[1:5,1:5]

### lets check everything is ok
dim(mrnaNorm)
dim(mrnaIDs)
mrnaIDs5 = mrnaIDs[,c(1,2,3,4,5)] # extract first 5 patientIDs
mrnaNorm5x5 = mrnaNorm[1:5, 1:5] # first 5 rows and columns 
head(mrnaIDs5, 2) # display first two rows
head(mrnaNorm, 2) # display first two rows
summary(mrnaNorm5x5) # summary statistics

"The following instruction analyzes the patient IDs because the IDs can tell us whether a 
particular column is from a normal sample or a tumor sample. TCGA has data both from tumor 
samples and from normal samples taken from the same patient. This information can be found 
in the label of the column. For example in ID: TCGA.3C.AAAU.01A.11R.A41B.07 the type of
sample is indicated in the 4th group: 01A. Tumor types range from 01 - 09, normal types 
from 10 - 19 and control samples from 20 - 29. Therefore we are extracting in the cell 
below this code indicating whether the sample is from normal or not. 
This group is the fourth in the ID character string, obtained with 'strsplit', 
and we get the two characters of interest by taking the first two characters with 'substr'.
We do this for all the IDs by applying with 'lapply' the same process to all the IDs
in 'mrnaIDs'."

samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
# extract the sample type (tumor / normal), excluding the first column since it is a gene ID
sampleType <- as.data.frame(samp)
dim(sampleType)
sampleType[1:5,1:5]

"Next we are going to count how many normals and how many tumor samples we have. 
The function 'unique' gives us the different types we have. 
The function 'table' calculates how many samples we have of each type. 
We see that we have 1100 of type '01' or '06' for tumor, and 112 of type '11' for normal."

unique(sampleType)
# extracts how many unique objects there are 

tab <- table(unlist(sampleType))
tab

"We are going to associate a class of '1' for the tumor samples and of '0' for the normal 
samples. We check that we have 112 normals (class '0') and 1100 tumors (class '1').
We do this for all the types by applying with 'lapply' the same process to all the types in
'samp'."


sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaClass)
table(unlist(sampClass))


# Asigning Numerical 0 and 1 class
sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum) 
mrnaClassNum

"Here I provide a function that calculates a ranking of features (genes in this case) 
with the BSS/WSS method. You do not need to understand all the details - only to 
execute the following cell. This will make the bssWssFast function available for
later processing."

X = t(mrnaNorm[, -1])
class(X)
y <- X[1:5,1:5]
class(y)
as.data.frame(lapply(y, as.numeric))

summary(mrnaNorm[,-1][1:5,1:5])

bssWssFast <- function (X, givenClassArr, numClass=2)
  # between squares / within square feature selection
  {
  X <-data.frame(lapply(X, as.numeric))
  classVec <- matrix(0, numClass, length(givenClassArr))
  for (k in 1:numClass) {
    temp <- rep(0, length(givenClassArr))
    temp[givenClassArr == (k - 1)] <- 1
    classVec[k, ] <- temp
  }
  classMeanArr <- rep(0, numClass)
  ratio <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    overallMean <- sum(X[, j]) / length(X[, j])
    for (k in 1:numClass) {
      classMeanArr[k] <- 
        sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
    }
    classMeanVec <- classMeanArr[givenClassArr + 1]
    bss <- sum((classMeanVec - overallMean)^2)
    wss <- sum((X[, j] - classMeanVec)^2)
    ratio[j] <- bss/wss
  }
  sort(ratio, decreasing = TRUE, index = TRUE)
}

"We invoke the bssWSSFast function by passing to it the dataframe of gene expressions
'mrnaNorm' (all numeric, after excluding the first column, which is not numeric since
it contains gene names) and another dataframe of classes associated with each sample. 
Note that we also transpose the 'mrnaNprm' dataframe because we watb to select genes, 
therefore the genes need to be feature / columns, while in 'mrnaNorm' genes are rows"


summary(mrnaNorm[2, ])
bss <- bssWssFast(t(mrnaNorm[, -1]), t(mrnaClassNum), 2)
# returns a list with two elements: a list of bss/wss ratios in decreasing order bss[1] and a list of features bss[2]

bss$ix[1:50]

"We can then list the genes by their name as provided in the 'mrnaNorm' 
dataset, where the first column is the list of gene names. The names have two 
codes separated by character '|':
HUGO code Entrez code
For example, FHL1|2273 is gene with HUGO code FHL1 and Entrez code 2273.
We will continue working on these genes in a future module."

genes <- mrnaNorm[bss$ix[1:50],1]
genes

##Heatmap
"Heatmaps are visualizations that represent the level of a variable across different 
samples or groups using different colors to represent levels.
In this particular example, we would like to represent the different expression levels of 
some of the genes selected above between normal and tumor samples.
We will represent tumor samples in red and normals as blue at the top.
First we are going to select the top 100 genes ranked, and associate the class values 0/1 
in a row at the bottom of mrnaSetClass.
Because selection by column is easier, we transpose the dataframe (inverse tows and columns)
We then separate the normals and the tumors in two datasets, selecting only 112 tumors
so that we have an equal number of normals and tumors.
We then combine the two sets together.
We draw the heatmap by selecting several options:
1) we choose a color scheme as red for tumors and blue for normals for the top bar.
2) we apply the color scheme to the classes in the 'both' dataframe.
3) we draw the heatmp with options 'scale = "row"' to normalize the rows (gene values), 'col=topo.colors(100)' as a color scheme for the heatmap cells, and 'ColSideColors=groupColors' to specify the colors chosen for the side bars.
What we find on the heatmap is that it finds several groupins in the patients, mostly two groups of normals, and two groups of tumors. So this is interesting. We can also see some gene expressions characterizing the groups."

# merge class values into mrnaNorm as the last row  - after the genes and select only top 100 genes
mrnaSetClass <- rbind(mrnaNorm[bss$ix[1:100],-1], setNames(mrnaClassNum, names(mrnaNorm[,-1])))
dim(mrnaSetClass)
# 101  1212       101 genes as rows and 1212 patients as columns
# to select by class, transpose the matrix
transpose <- as.data.frame(t(mrnaSetClass))
colnames(transpose)[101] <- "class"
# select normals
normals <- subset(transpose, transpose$class == 0)
dim(normals)
#[1] 112 101
# select tumors
tumors <- subset(transpose, transpose$class == 1)
dim(tumors)
#[1] 1100 101
# combine them together
both <- rbind(normals, tumors[1:112,])
dim(both)
# select a color scheme - red for tumor and blue for normal
color.map <- function(class) { if (class==0) "#FF0000" else "#0000FF" } # red for tumor and blue for normal
# create heatmap
groupColors <- unlist(lapply(both$class, color.map))
heatmap(as.matrix(t(both)), scale = "row", col=topo.colors(100), ColSideColors=groupColors)



