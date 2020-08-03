library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(ggplot2)

setwd("~/work/Bioinformatics/DataSets/")
dataset <- read.table("heart-ch.txt",
                      header = TRUE, sep=",", quote="")
summary(dataset)
dim(dataset)
nrow(dataset)
ncol(dataset)
dim(dataset)
head(dataset,0)  


#md.pattern(dataset)  # (mice package) displays all the missing values, NA for missing values
md.pattern(dataset)
# where '0' represents a missing value, and '1' represents a present value. Each row represents an increasing number of missing values, starting with 0 missing values, then 1, etc. (right column).
# The left column indicates how many rows have this number of missing values. Here for example, 297 rows have no missing values.
# From the 'VIM' package, we can display missing data in a different way with the 'aggr' function. 
# This function ranks the variables by decreasing number of missing values
mice_plot <- aggr(dataset, col=c("green","red"),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(dataset), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern")) # (VIM) display graphically missing values

mice_plot <- aggr(dataset, col=c("green","red"))

"To deal with missing values, one main method is to impute 
these values, for example with the mean of the column 
where the missing value is located. This can be accomplished for example with the 'impute' function. In the example below, 
we are imputing 'chol' with its mean, and copy the results in
a new variable called 'chol_imputed'"

dataset[,"chol_imputed"] <- with(dataset, impute(chol, mean))  # (Hmisc) impute with mean
summary(dataset)
  
"A function 'Normalize' exists in 'RWeka' package to normalize 
using ZScore normalization."

#dataset_n <- Normalize(chol_imputed ~., data = dataset) # (RWeka) normalizes all numeric except chol_imputed
#dataset_n <- Normalize(~ chol_imputed , data = dataset) # normalizes only chol_imputed
dataset_n <- Normalize(~. , data = dataset) # normalizes all variables
summary(dataset_n)
head(dataset_n)
## Descritization

dataset$chol_bin <- as.numeric(cut2(dataset$chol_imputed, g=3)) # (Hmisc) create 3 bins g = 3 quantile bins equal-depth
summary(dataset)
head(dataset)

"For an equal-width discretization, 5 bins can be created with
'cut' for example. We see in 'summary' that 'chol_bin' has five
values 1, 2, 3, 4, and 5, each representing a bin."
dataset[,"chol_bin"] <- as.numeric(cut(dataset[,"chol_imputed"], 5)) # create bins of same width
summary(dataset)

#Data Visualization##
"Data understanding is a data analytics step preceding data preprocessing.
However it can also be used at any time during analysis to better
understand the data. 
Function 'ggplot' below creates a scatterplot of the data 
in 'dataset' representing two numeric variables:
'chol' as a function of 'age' and using two different colors for 'num'. 
The plot shows more narrowing of arteries in older ages."

ggplot(dataset,aes(x=age,y=chol, color=num)) + geom_point(size = 4) # scatterplot (age, chol)

"Another plot can be a stacked histogram for a 
categorical variable ('chest_pain'). 
The following plot shows that some types of chest pain 
(asymptomatic pain) are more associated with heart disease
than others."
ggplot(dataset, aes(chest_pain, fill=factor(num))) + geom_bar() # stacked histogram for categorical variable

"A stacked histogram can also be used for numeric variables such as 'age'.
The following plot shows that heart disease is more prevalent when age increases."
ggplot(dataset, aes(x=age, fill=factor(num))) + geom_bar() # stacked histogram for numeric variable








