#!/usr/bin/env Rscript

#===============================================================================================
Usage<-function(){
cat("Usage: Rscript findDiffpeaks_hs.r [readcounts] [patientinfo] [adj_Pvalue] [fold_change]\n\n",

    "Parameters:\n",
    "[readcounts]     hsATACseq_readcounts.txt, including all 5 samples\n",
    "[patientinfo]    patientinfo.txt, containing information of all 5 samples\n",
    "[adj_Pvalue]     adjusted Pvalue for defining DARs, i.e. 0.05\n",
    "[fold_change]    fold change for defining DARs, i.e. 2\n\n",

    "Example: (R-3.6.0)\n",
    "Rscript findDiffpeaks_hs.r hsATACseq_readcounts.txt patientinfo.txt 0.05 2\n\n",

    "Function: Use DESeq2 to find differentially accessible regions (DARs) and generate plots. \"hs\" in findDiffpeaks_hs.r means human version.\n",
    "Qirui Zhang (qirui.zhang@med.lu.se)\n",
    "05-07-2021\n\n"
    )
}

args<-commandArgs(TRUE)
if(length(args)!=4){Usage();quit();}


cat("\n=========================================================================================\n")
# Load libraries
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Loading libraries and arguments...\n\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
pdf("diffPeaks_plots.pdf", useDingbats=FALSE)
options(scipen=20)

reads.count.file<-args[1]  # hsATACseq_readcounts.txt
patient.info.file<-args[2]  # patientinfo.txt
adjust.pvalue<-as.numeric(args[3])  # 0.05
fold.change<-as.numeric(args[4])  # 2

# read data
cat("reading files...\n\n")
reads.count<-read.table(reads.count.file, header=T, stringsAsFactors=F)
rownames(reads.count)<-reads.count$peakID
reads.count<-reads.count[,5:9]
colnames(reads.count)<-c("1D", "3D", "46D", "48D", "50D")

patient.info<-read.table(patient.info.file, header=T, stringsAsFactors=F)
patient.info$Month<-as.factor(patient.info$Month)
patient.info$FusionGene<-as.factor(patient.info$FusionGene)
patient.info$ActivatingMut<-as.factor(patient.info$ActivatingMut)


cat("\n=========================================================================================\n")
# Run DESeq2
cat("Running differential comparison...\n\n")

# compare month "morethan3" to "lessthan3"
cat("Comparing morethan3 to lessthan3 month...\n")
dds.month<-DESeqDataSetFromMatrix(countData=reads.count, colData=patient.info, design= ~ Month)
dds.month<-dds.month[rowSums(counts(dds.month))>=10,]
dds.month$Month<-relevel(dds.month$Month, ref="lessthan3")
dds.month<-DESeq(dds.month)
res.month<-results(dds.month, contrast=c("Month", "morethan3", "lessthan3"), alpha=adjust.pvalue)
summary(res.month)

up_fc2<-res.month[!is.na(res.month$padj) & res.month$padj<adjust.pvalue & res.month$log2FoldChange>=log2(fold.change),]
up_fc2.num<-nrow(up_fc2)
down_fc2<-res.month[!is.na(res.month$padj) & res.month$padj<adjust.pvalue & res.month$log2FoldChange<= -log2(fold.change),]
down_fc2.num<-nrow(down_fc2)
dar.peakid.month<-c(rownames(up_fc2), rownames(down_fc2))
cat("Increased peaks with |foldchange|>2:", up_fc2.num, "\n")
cat("Decreased peaks with |foldchange|>2:", down_fc2.num, "\n\n")
write.table(as.data.frame(up_fc2), "diffPeaks_MorevsLess3month_up.tsv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(down_fc2), "diffPeaks_MorevsLess3month_down.tsv", row.names=T, col.names=T, quote=F, sep="\t")

# compare activating mutations to non-activating mutations
cat("Comparing activating mutations to non-activating mutations...\n")
dds.activating<-DESeqDataSetFromMatrix(countData=reads.count, colData=patient.info, design= ~ ActivatingMut)
dds.activating<-dds.activating[rowSums(counts(dds.activating))>=10,]
dds.activating$ActivatingMut<-relevel(dds.activating$ActivatingMut, ref="NONE")
dds.activating<-DESeq(dds.activating)
res.activating<-results(dds.activating, contrast=c("ActivatingMut", "RAS", "NONE"), alpha=adjust.pvalue)
summary(res.activating)

up_fc2<-res.activating[!is.na(res.activating$padj) & res.activating$padj<adjust.pvalue & res.activating$log2FoldChange>=log2(fold.change),]
up_fc2.num<-nrow(up_fc2)
down_fc2<-res.activating[!is.na(res.activating$padj) & res.activating$padj<adjust.pvalue & res.activating$log2FoldChange<= -log2(fold.change),]
down_fc2.num<-nrow(down_fc2)
dar.peakid.activating<-c(rownames(up_fc2), rownames(down_fc2))
cat("Increased peaks with |foldchange|>2:", up_fc2.num, "\n")
cat("Decreased peaks with |foldchange|>2:", down_fc2.num, "\n\n")
write.table(as.data.frame(up_fc2), "diffPeaks_ActivatingvsNone_up.tsv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(down_fc2), "diffPeaks_ActivatingvsNone_down.tsv", row.names=T, col.names=T, quote=F, sep="\t")

# normalize reads count and vst transform
cat("Normalizing reads count ...\n")
normalized.counts=as.data.frame(counts(dds.activating, normalized=TRUE))
normalized.counts$mean=apply(normalized.counts, 1, mean)
write.table(as.data.frame(normalized.counts), "atac-seq_deseq2norm_hs.tsv", row.names=T, col.names=T, quote=F, sep="\t")
vst.month=varianceStabilizingTransformation(dds.month, blind=FALSE)
vst.activating=varianceStabilizingTransformation(dds.activating, blind=FALSE)
write.table(as.data.frame(assay(vst.activating)), "atac-seq_vst_hs.tsv", quote=F, row.names=T, col.names=T, sep="\t")


cat("\n=========================================================================================\n")
cat("Correlation heatmap of samples...", "\n")
sampleDist<-dist(t(assay(vst.activating)))
sampleDist.mx<-as.matrix(sampleDist)
rownames(sampleDist.mx)<-patient.info$Patient
colnames(sampleDist.mx)<-patient.info$Patient
colors<-rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100))
pheatmap(as.data.frame(sampleDist.mx), clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist, color=colors, show_rownames=T)

### Not Use (below) ###########################################################
if(FALSE){
# Plot PCA using DESeq2::plotPCA
cat("PCA plot using DESeq2::plotPCA ...", "\n")
PCAbyDESeq2<-function(vst, group){
    data<-plotPCA(vst, intgroup=group, returnData=TRUE)
    percentVar<-round(100*attr(data, "percentVar"))
    pca.plot<-ggplot(data, aes(PC1, PC2, color=group))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_rect(size=1),text=element_text(size=18),axis.text=element_text(size=15),legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=1)
    out1=pca.plot
    out2=pca.plot+geom_text_repel(aes(label=rownames(data)), size=3)
    out=list(out1, out2)
    return(out)
}

PCAbyDESeq2(vst.month, "Month")
PCAbyDESeq2(vst.activating, "ActivatingMut")
}
### Not Use (above) #########################################################


# Plot PCA using in-house script (recommended!)
cat("PCA plot using in-house script ...", "\n")
PCAinHouse<-function(vst, group){
    vst=as.data.frame(assay(vst))
    vst=t(vst)
    pca.results=prcomp(vst, scale=TRUE)
    percentVar=round(((pca.results$sdev^2)/sum(pca.results$sdev^2))*100, 1)
	scree.plot=barplot(percentVar, main=paste("Scree plot of ", group, sep=""), xlab="Principal component", ylab="Variation percent")
    data=as.data.frame(pca.results$x)
    data$group=patient.info[,group]
    pca.plot=ggplot(data, aes(PC1, PC2, color=group))+geom_point(size=4)+scale_color_manual(name="Group", values=c("#F8766D", "#00BFC4"))+xlab(paste("PC1: ", percentVar[1], "% variance", sep=""))+ylab(paste("PC2: ", percentVar[2], "% variance", sep=""))+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(color="black",size=0.5), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12))+coord_fixed(ratio=1)
    out1=pca.plot
    out2=pca.plot+geom_text_repel(aes(label=rownames(data)), size=3)
    out=list(out1, out2)
    return(out)
}

PCAinHouse(vst.month, "Month")
PCAinHouse(vst.activating, "ActivatingMut")


cat("\n=========================================================================================\n")
cat("Volcano plots ...", "\n\n")
VolcanoPlot<-function(res, plot.title){
    volcano<-as.data.frame(res)
    volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange)>=1, ifelse(volcano$log2FoldChange>=1, "Up", "Down"), "No"))
    volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)
    volcano$group<-rep(0)
    volcano[which(volcano$padj >= adjust.pvalue | (volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange) < 1)),8]=rep(1)
    volcano[which(volcano$padj < adjust.pvalue & volcano$log2FoldChange < -1),8]=rep(2)
    volcano[which(volcano$padj < adjust.pvalue & volcano$log2FoldChange > 1),8]=rep(3)

    p<-ggplot(volcano, aes(log2FoldChange, -log10(padj)))+geom_point(data=volcano[which(volcano$group==1),], color="gray55", alpha=0.7)+geom_point(data=volcano[which(volcano$group==2),], color="royalblue3", alpha=0.75)+geom_point(data=volcano[which(volcano$group==3),], color="firebrick3", alpha=0.75)
    p+labs(title={plot.title}, x="log2FoldChange", y="-log10(padj)")+geom_hline(yintercept=-log10(adjust.pvalue), linetype=2,color="black")+geom_vline(xintercept=c(-log2(2), log2(2)), linetype=2, color="black")+theme_bw()+theme(text=element_text(size=12), panel.grid=element_blank(), panel.border=element_rect(size=0.5), axis.line=element_line(color="black", size=0.8), axis.ticks=element_line(color="black", size=0.5), axis.text=element_text(color="black"), legend.position = c(.95, .95))
}

VolcanoPlot(res.month, "Morethan3 vs Lessthan3")
VolcanoPlot(res.activating, "ActivatingMut vs None")


cat("\n=========================================================================================\n")
cat("Volcano plots in modified version ...", "\n\n")
VolcanoPlotV2<-function(res, plot.title){
    volcano<-as.data.frame(res)
    volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange)>=1, "yes", "no"))
    volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)

    p<-ggplot(volcano, aes(x=log2FoldChange, y=log2(baseMean)))+geom_point(data=subset(volcano, significant=="no"), color="gray55")+geom_point(data=subset(volcano, significant=="yes"), aes(color=padj), alpha=0.8)+scale_color_gradientn("Padj", colors=c(rev(brewer.pal(9, "OrRd"))[c(2,5,6,7,8,9)]), limits=c(0, adjust.pvalue), breaks=c(seq(0.01, adjust.pvalue, 0.02)), labels=as.character(c(seq(0.01, adjust.pvalue, 0.02))))
    p+labs(title={plot.title}, x="Log2foldchange", y="Log2basemean")+geom_vline(xintercept=c(-1,1), linetype=2,color="black")+theme_bw()+theme(text=element_text(size=12), axis.title=element_text(size=12), axis.text=element_text(color="black"), panel.grid=element_blank(), panel.border=element_blank(), axis.line=element_line(color="black", size=0.8), axis.ticks=element_line(color="black", size=0.5), legend.text=element_text(size=12), legend.position=c(0.9, 0.85))+coord_flip(xlim=c(-5, 5))
}

VolcanoPlotV2(res.month, "Morethan3 vs Lessthan3")
VolcanoPlotV2(res.activating, "ActivatingMut vs None")


cat("\n=========================================================================================\n")
cat("Heatmap of DARs ...", "\n\n")
# morethan3 vs lessthan3
vst.month.dar<-vst.month[dar.peakid.month,]
vst.month.dar<-as.data.frame(assay(vst.month.dar))
anno.label<-data.frame(ActivatingMut=c("lessthan3", "morethan3", "morethan3", "lessthan3", "lessthan3"))
rownames(anno.label)<-patient.info$Patient
anno.color<-list(ActivatingMut=c("lessthan3"="#F8766D", "morethan3"="#00BFC4"))
pheatmap(vst.month.dar, main="morethan3 vs lessthan3", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)

# activatingMut vs none
vst.activating.dar<-vst.activating[dar.peakid.activating,]
vst.activating.dar<-as.data.frame(assay(vst.activating.dar))
anno.label<-data.frame(ActivatingMut=c("RAS", "RAS", "NONE", "RAS", "NONE"))
rownames(anno.label)<-patient.info$Patient
anno.color<-list(ActivatingMut=c("RAS"="#F8766D", "NONE"="#00BFC4"))
pheatmap(vst.activating.dar, main="RAS vs NONE", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)

dev.off()


#===============================================================================================
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Done with analysis!", "\n\n")

