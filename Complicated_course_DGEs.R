#Author : Shayantan Banerjee
#assuming the clinical_and_genomic.txt file is in your path and loaded into the R workspace with the filename data
arr=numeric(276)
arr[which(data$Complicated.Course=="YES")]="YES"
arr[which(data$Complicated.Course=="NO")]="NO"
#tail(clinical_data_without_controls)
arr[which(arr=="0")]<-"X"
arr[which(arr=="YES")]<-"1"
arr[which(arr=="NO")]<-"0"
gsms=paste(arr, collapse = "")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sel <- which(sml != "X")
sml <- sml[sel]
x=newrma_aggregrated[,-c(1)][,sel]
dim(x)
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
g=t(x)
g$description <- fl
arr=arr[which(arr=="0" | arr=="1")]
arr
design <- model.matrix(~0+factor(c(arr)))
colnames(design) <- levels(fl)
#x1=newrma_aggregrated[,-c(1)][,sel]
x[] <- lapply(x, function(x) as.numeric(as.character(x)))
#str(x1)
x1=as.data.frame(x)
fit <- lmFit(x, design)
contrast.matrix <- makeContrasts(diff=G1-G0, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tT= topTable(fit2, number=Inf, p.value = 0.05, adjust.method = "fdr", sort.by="logFC")
head(tT)
tT=cbind(tT,Gene_name=newSymbols_new[match(rownames(tT),newSymbols_new$PROBEID),4])
#tT=tT[!duplicated(tT$ID),]
#sorting by adjusted p value
tT_pvalue=tT[order(tT$adj.P.Val,decreasing = FALSE),]
tT_pvalue_0.05=tT_pvalue[which(tT_pvalue$adj.P.Val<0.05),] #756 genes
#considering a log fold change of > 1
tT_signi_up_down=subset(tT_pvalue, tT_pvalue$logFC > 1 | tT_pvalue$logFC < -1) #18 genes
#plotting heatmaps
heatmap.plus(as.matrix(exp1),labCol = FALSE,scale = "row",col = bluered(20),margin=c(20,13),ColSideColors = myCols,,trace="none", density="none")legend(0.5,0.1,legend=c("Complicated","Non Complicated"),fill=c("blue","yellow"),cex=0.5)
legend(0.1,0.1,legend=c("Survivor","Non survivor"),fill=c("black","purple"),cex=0.5)
#plotting volcanoplots
EnhancedVolcano(tT_signi_up_down, lab = tT_signi_up_down$ID,
                x = 'logFC',
                y = 'adj.P.Val',pCutoff = 0.05,FCcutoff = 0.5
              ,transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                title= "Complicated vs uncomplicated without controls")
