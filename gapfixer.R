#!/usr/bin/env Rscript
## Usage is gapfixer.R <ref_file> <medaka_file> <final_file> 

args = commandArgs(trailingOnly=TRUE)

library(DECIPHER)

ref_file = args[1]
medaka_file = args[2]
final_file = args[3]

ref_fasta = readDNAStringSet(ref_file)
medaka_fasta = readDNAStringSet(medaka_file)

final_fasta = DNAStringSet()

for(i in seq(1, length(ref_fasta), 1)){
    #print(i)
    ref_seq = ref_fasta[i]
    medaka_pos = pmatch(names(ref_seq)[1], names(medaka_fasta))
    if(!is.na(medaka_pos)){
        #print("Yeah")
        medaka_seq = medaka_fasta[medaka_pos]
        align_seq = c(ref_seq, medaka_seq)
        aligned <- AlignSeqs(align_seq)
        #BrowseSeqs(aligned)
       
       
        ref_letter <- unlist(strsplit(as.character(aligned[[1]]), split = ""))
        medaka_letter <- unlist(strsplit(as.character(aligned[[2]]), split = ""))
        #Build based on Reference Gap
        #snps <- ((ref_letter != medaka_letter) & (medaka_letter != "-") & (ref_letter != "-"))
        #ref_letter[snps] <- medaka_letter[snps]
        #ref_letter[ref_letter == "-"] <- ""
        #out_seq <- paste(ref_letter, collapse = "")
       
        #Build based on Medaka Gap
        medaka_letter[ref_letter == "-"] <- ""
	medaka_letter[medaka_letter == "-"] <- ref_letter[medaka_letter == "-"]
        out_seq <- paste(medaka_letter, collapse = "")
        final_sequence <- DNAStringSet(out_seq)
        names(final_sequence) <- names(ref_seq)
        final_fasta <- c(final_fasta, final_sequence)
    }
   
    writeXStringSet(final_fasta, final_file )    
}
