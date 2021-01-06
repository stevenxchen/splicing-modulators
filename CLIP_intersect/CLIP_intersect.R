# Intersect CLIPseq RBP-binding data with splicing annotations ####

rm(list=ls())
setwd("/Data")

# Make BED file for each splice type GTF File ####
# Only need to do this once for each of SE/A5SS/A3SS/MXE/RI
### Read in splice events and make BED file ###

splice.type <- "RI"
splice_events <- read.table(file=paste0("../gencode.v19.",splice.type,".txt"), header = TRUE, sep = "", stringsAsFactors = FALSE)

rn <- nrow(splice_events)
splice_bed <- data.frame(chrom=character(rn), chromStart=numeric(rn), chromEnd=numeric(rn), name=character(rn),
                         geneSymbol=numeric(rn), strand=character(rn),
                         stringsAsFactors = FALSE)

# Populate the table
splice_bed$chrom <- splice_events$chr

if (splice.type == "SE"){
  splice_bed$chromStart <- splice_events$exonStart_0base
  splice_bed$chromEnd <- splice_events$exonEnd
} else if (splice.type == "RI"){
  splice_bed$chromStart <- splice_events$riExonStart_0base
  splice_bed$chromEnd <- splice_events$riExonEnd
} else if (splice.type == "A5SS" || splice.type == "A3SS"){
  splice_bed$chromStart <- splice_events$longExonStart_0base
  splice_bed$chromEnd <- splice_events$longExonEnd
} else if (splice.type == "MXE"){
  splice_bed$chromStart <- splice_events$X1stExonStart_0base
  splice_bed$chromEnd <- splice_events$X2ndExonEnd
} else{
    print("Splice type incorrect")
}

splice_bed$name <- splice_events$GeneID
splice_bed$strand <- splice_events$strand
splice_bed$geneSymbol <- splice_events$geneSymbol
splice_bed$ID <- splice_events$ID

# Add 300 upstream and 300 downstream base positions
splice_export <- splice_bed
splice_export$chromStart <- splice_export$chromStart -300
splice_export$chromEnd <- splice_export$chromEnd +300

# Export file and run bedtools
write.table(splice_export, file=paste0(splice.type,"-splice-sites-300.bed"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

# Output:
# SE-plice-sites-300.bed
# A5SS-splice-sites-300.bed
# A3SS-splice-sites-300.bed
# MXE-splice-sites-300.bed
# RI-splice-sites-300.bed


# SE Bedtools Intersect with CLIP seq binding sites ####

# bedtools intersect -a SE-splice-sites-300.bed -b SRSF1* -u -s > SRSF1-SE-intersect
# -u for unique
# -s for strandedness

# # Load GTF file for SE events
# splice_events <- read.table(file="../fromGTF.SE.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

rbp <- "SRSF1" # Define RBP

# Load event reads from dataset ####
load("SE.event_reads.Rdata")

outDir <- "/Output/" # Define output directoru

bedArg <- paste0("bedtools intersect -a SE-splice-sites-300.bed -b ", CLIPdir, "ENCFF530MIN.bed ",
                 CLIPdir, "ENCFF788UFQ.bed ", "-u -s > ", outDir, rbp[i], "-SE.bed")
system(bedArg)

setwd(outDir)
file <- paste0(rbp[i], "-SE.bed")
fileIn <- read.delim(file = file, header = FALSE, stringsAsFactors = FALSE)

# set column names
colnames(fileIn) <- c("chrom", "chromStart", "chromEnd", "name", "geneSymbol", "strand", "ID")

# Remove 300bp upstream and downstream
fileIn$chromStart <- fileIn$chromStart +300
fileIn$chromEnd <- fileIn$chromEnd -300

outReads <- event_reads[fileIn$ID,]
outName <- fileIn[c("name", "geneSymbol")]
                  
save(outReads, outName, file=paste0(outDir, rbp[i], ".SE-reads.Rdata"))
  