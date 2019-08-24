#Author : Akram Mohammed and Rishikesan Kamaleswaran and Shayantan Banerjee
install.packages("magrittr")
install.packages("dplyr")
install.packages("AnnotatedDataFrame")
#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(BiocUpgrade)
# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hgu133plus2cdf")
biocLite("hgu133plus2probe")
biocLite("hgu133plus2.db")
library(limma)
library(magrittr)
#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(hgu133plus2.db)
library(dplyr)
library(AnnotationDbi)
#Set working directory for download
#Download GSE66099 raw data from GEO and place it in your current working 
setwd("/home/shayantan/")
untar("GSE66099_RAW.tar",exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("/home/shayantan/")
pd<-read.AnnotatedDataFrame("Complicated_course_label.csv",sep=",", header=T)
setwd("/home/shayantan/data")
rawData=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hgu133plus2", phenoData=pd) #From bioconductor
class(normData)

#perform RMA normalization
normData=rma(rawData)
#filterig expressin set data
#ans <- featureFilter(normData)
#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(normData)
dim(rma) #10093 by 276
#Filter Affy controls
rma=rma[which(!grepl("AFFX", rownames(rma))),]

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe set identifiers to Entrez gene symbols and IDs and then combine with raw data.
#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hgu133plus2.db") #Annotations at the exon probeset level

probes=rownames(rma)


collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
  
}

# Redefinition of ae.annots
newSymbols <- AnnotationDbi::select(
  x       = hgu133plus2.db,
  keys    = rownames(rma),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarise_all(funs(collapser)) %>%
  ungroup

fd <- new("AnnotatedDataFrame",
          data = data.frame(rma, stringsAsFactors = FALSE)
)

#The next step requires you to have the clinical data file with the sample names and the complicated course labels without control samples 
rma_new=rma[,c(which(colnames(rma) %in% clinical_data_without_controls$Samples))]
rownames(fd) <- newSymbols$PROBEID
newrma=cbind(newSymbols,rma_new)
#removing probes that didn't match any known genes
newrma=newrma[-c(which(newrma$SYMBOL=="")),]


#preparing  the final expression matrix with the gene symbols as the row names
#newrma=cbind(newSymbols_new,rma_new)
newrma_numeric <- lapply(newrma[,-c(1,2,3,4)], function(x) as.numeric(as.character(x)))
newrma_numeric=as.data.frame(newrma_numeric)
newrma_numeric['SYMBOL']=newrma$SYMBOL
rownames(newrma_numeric)=newrma$PROBEID
#In case of multiple probes match to the same gene id, the expression values are averaged
newrma_aggregrated=aggregate(newrma_numeric[, -c(229)],
          by = list(Gene = newrma$SYMBOL),
          FUN = mean,
          na.rm = TRUE)
probeid=as.data.frame(newSymbols_new[match(newrma_aggregrated$Gene,newSymbols_new$SYMBOL),1])
rownames(newrma_aggregrated)=probeid$PROBEID
dim(newrma_aggregrated) #21724 by 229
