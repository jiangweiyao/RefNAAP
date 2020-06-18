#!/usr/bin/env Rscript
## Usage is scaffold_cutter.R <topN> <min_cov> <fasta_file> <depth_file> <out_fasta> <out_table>

library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

topN = as.numeric(args[1])
min_cov = as.numeric(args[2])
fasta_file = args[3]
depth_file = args[4]
out_fasta = args[5]
out_table = args[6]


depth_count <- read.table(file = depth_file, sep = "\t", stringsAsFactors = FALSE)
colnames(depth_count) <- c("RefSeq", "Position", "Depth")


depth_count <- subset(depth_count, Depth > 0)

#Find beginning boundary and mark as false
depth_count$Begin <- (depth_count$Position -1 == head(c(-1, depth_count$Position), -1)) & (depth_count$RefSeq == head(c("-1", depth_count$RefSeq), -1))
depth_count$End <- (depth_count$Position + 1 == c(depth_count$Position[-1], -1)) & (depth_count$RefSeq == c(depth_count$RefSeq[-1], "-1"))

depth_Begin <- subset(depth_count, (depth_count$Begin == FALSE ))
depth_End <- subset(depth_count, (depth_count$End == FALSE ))

range_list  <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(range_list) <- c("RefSeq", "Begin", "End", "AveDepth")

for(i in seq(1, dim(depth_Begin)[1])){
    RefSeqi=as.character(depth_End[i,"RefSeq"])
    RefSeqii=as.character(depth_Begin[i,"RefSeq"])
    #print(RefSeqi == RefSeqii)
    Begini=as.numeric(depth_Begin[i,"Position"])
    Endi=as.numeric(depth_End[i,"Position"])
    AveDepthi=mean(subset(depth_count, (RefSeq == RefSeqi)&(Position >= Begini)& (Position <= Endi))$Depth)
    range_summary <- data.frame(RefSeq=RefSeqi,
                                Begin=Begini,
                                End=Endi,
                                AveDepth=AveDepthi)
    range_list <- rbind(range_list, range_summary)
}

write.csv(range_list, paste0(out_table, ".raw"), row.names = FALSE)

range_list <- subset(range_list, AveDepth > min_cov)

RefSeq_coverage <- aggregate(range_list$AveDepth, by=list(range_list$RefSeq), FUN=mean)
RefSeq_coverage <- RefSeq_coverage[order(order(RefSeq_coverage[,2], decreasing = TRUE)),]

topN_ref <- as.character(head(RefSeq_coverage, topN)[,1])
range_list_topN <- subset(range_list, RefSeq %in% topN_ref)
write.csv(range_list_topN, out_table, row.names = FALSE)
   
fasta <- readDNAStringSet(fasta_file)
fasta_list <- DNAStringSet()

for(i in seq(1, dim(range_list_topN)[1], 1)){
    fastai <- DNAStringSet(fasta[[as.character(range_list_topN[i, 1])]][range_list_topN[i, 2]:range_list_topN[i, 3]])
    #print(i)
    fasta_list <- c(fasta_list, fastai)
    #print(fasta_list)
    #writeXStringSet(fastai, "fasta.fasta", append=TRUE)
}

names(fasta_list) <- paste(as.character(range_list_topN[, 1]), as.character(range_list_topN[, 2]), as.character(range_list_topN[, 3]), sep = "_")
writeXStringSet(fasta_list, out_fasta)
