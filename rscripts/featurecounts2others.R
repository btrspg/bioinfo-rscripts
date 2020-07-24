#! /usr/bin/env Rscript

library("DESeq2")

args<-commandArgs(TRUE)
if(length(args) != 3){
    stop("featurecounts2others.R <> <clinical_info> <prefix> ")
}



featurecount_list <- args[1]
clinical_info <- args[2]
prefix <- args[3]

depth.filter <- 10
qvalue.filter <- 0.1
pvalue.filter <- 0.1
ct <- read.table(count_table,sep="\t",header=TRUE,row.names=1)
ci <- read.table(clinical_info,sep="\t",header=TRUE,row.names=1)
dds <- DESeqDataSetFromMatrix(countData = ct,
                              colData = ci,
                              design = ~ batch + condition)
keep <- rowSums(counts(dds)) >= depth.filter
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = control)
dds <- DESeq(dds)
res <- results(dds)
cat(summary(res))
dir.create(dirname(prefix))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="{prefix}_ordered.csv")
resSig_q <- subset(resOrdered, padj < qvalue.filter)
resSig_p <- subset(resOrdered, pvalue < pvalue.filter)
rld <- rlog(dds, blind=FALSE)
write.csv(assay(rld),file="{prefix}_rld.csv")