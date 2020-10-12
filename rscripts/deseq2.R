#! /usr/bin/env Rscript

library("DESeq2")

args<-commandArgs(TRUE)
if(length(args) != 4){
    stop("deseq2.R <featurecount> <clinical_info> <control> <prefix> ")
}



count_table <- args[1]
clinical_info <- args[2]
control <- args[3]
prefix <- args[4]

depth.filter <- 10
qvalue.filter <- 0.1
pvalue.filter <- 0.05
ct <- read.table(count_table,sep="\t",header=TRUE,row.names=1,check.names=F)
ci <- read.table(clinical_info,sep="\t",header=TRUE,row.names=1,check.names=F)
dds <- DESeqDataSetFromMatrix(countData = ct,
                              colData = ci,
                              design = ~ batch + condition)
keep <- rowSums(counts(dds)) >= depth.filter
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = control)
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=paste(prefix,"_ordered.csv",sep=''),quote = FALSE)
resSig_q <- subset(resOrdered, padj < qvalue.filter)
resSig_p <- subset(resOrdered, pvalue < pvalue.filter)
rld <- rlog(dds, blind=FALSE)
write.csv(assay(rld),file=paste(prefix,"_rld.csv",sep=''),quote = FALSE)