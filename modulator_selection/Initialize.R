# Regulator Initialize ####
rm(list=ls())
setwd(".Data/")
library(dplyr)

# Load Modulators ####
load("ROSMAP.imputed.expression.Rdata")

# Apply filtering criteria to modulators
rowSD <- apply(imputed.genes, 1, sd)
modulator <- modulator[rowSD > 0.1,]
modulator.gene.info <- gene_info[rowSD > 0.1,]

# Discretize modulator values into tertiles ####

Modulator.discrete <- as.data.frame(modulator)
samples <- ncol(Modulator.discrete)
for(i in 1:nrow(Modulator.discrete)){
  mod_temp <- as.numeric(Modulator.discrete[i,]) #Modulator Expression
  
  ### Top 1/3 values transformed to 1 and bottom 1/3 values transformed to 0 ###
  mod_rank <- rank(mod_temp)
  mod_low <- which(mod_rank <= samples*1/3)
  mod_mid <- which(mod_rank <= samples*2/3 & mod_rank > samples*1/3)
  mod_up <- which(mod_rank > samples*2/3)
  
  mod_discrete <- mod_temp
  mod_discrete[mod_up] <- 1
  mod_discrete[mod_mid] <- NA
  mod_discrete[mod_low] <- 0
  
  Modulator.discrete[i,] <- mod_discrete
}
save(Modulator.discrete, file="../Input/Modulator.discrete.Rdata")

### Discretize splicing counts into tertiles ####

rm(list=ls())
setwd("/Data/")
rbp <- "SRSF1" # define RBP

# SE event ####
load(file=paste0(rbp, ".SE-reads.Rdata"))
len <- ncol(outReads)

# Filter samples by requiring at least 10 reads of inclusion or exclusion 
# Filter that an event should be supported by at least 200 samples
# Filter splice targets with IQR(psi) > 0.1
incl_reads <- outReads[, seq(1,len,2)]
excl_reads <- outReads[, seq(2,len,2)]
total_reads <- incl_reads + excl_reads
denom <- incl_reads + 2*excl_reads
psi <- incl_reads / denom

excl_support <- rowSums(excl_reads >= 1) >= 1
table(excl_support)

total_support <- rowSums(total_reads >= 10) >= 1
table(total_support)

keep <- (excl_support + total_support) == 2
table(keep)

psi <- psi[keep,]
gene_id <- outName[keep,] # Make sure gene_id / outName stays consistent

iqr <- apply(psi, 1, function(x) IQR(x, na.rm=TRUE))
psi <- psi[iqr>0.10, ] # IQR of PSI > 0.10
gene_id <- gene_id[iqr>0.10, ]

num_supporting_samples <- numeric()

# save(psi, gene_id, file="SE.psi.Rdata")

for(i in 1:nrow(psi)){
  psi_rank <- as.numeric(rank(psi[i,], na.last="keep"))
  length_samples <- sum(!is.na(psi_rank))
  psi_low_tert <- which(psi_rank <= length_samples*1/3)
  psi_mid_tert <- which(psi_rank <= length_samples*2/3 & psi_rank > length_samples*1/3)
  psi_up_tert <- which(psi_rank > length_samples*2/3)
  
  num_supporting_samples <- c(num_supporting_samples, length_samples)
  
  psi[i, psi_up_tert] <- 1
  psi[i, psi_mid_tert] <- NA
  psi[i, psi_low_tert] <- 0
}

save(psi, gene_id, num_supporting_samples, file=paste0("../Input/", rbp, ".SE.PSI.discrete.Rdata"))

# A5SS event ####
load(file=paste0(rbp, ".A5SS-reads.Rdata"))

# Filter samples by requiring at least 10 reads of inclusion or exclusion 
# Filter that an event should be supported by at least 200 samples
# Filter splice targets with IQR(psi) > 0.1
incl_reads <- outReads[, seq(1,len,2)]
excl_reads <- outReads[, seq(2,len,2)]
total_reads <- incl_reads + excl_reads
denom <- incl_reads + 2*excl_reads
psi <- incl_reads / denom

excl_support <- rowSums(excl_reads >= 1) >= 1
table(excl_support)

total_support <- rowSums(total_reads >= 10) >= 1
table(total_support)

keep <- (excl_support + total_support) == 2
table(keep)

psi <- psi[keep,]
gene_id <- outName[keep,] # Make sure gene_id / outName stays consistent

iqr <- apply(psi, 1, function(x) IQR(x, na.rm=TRUE))
psi <- psi[iqr>0.10, ] # IQR of PSI > 0.10
gene_id <- gene_id[iqr>0.10, ]

num_supporting_samples <- numeric()

# save(psi, gene_id, file="A5SS.psi.Rdata")

for(i in 1:nrow(psi)){
  psi_rank <- as.numeric(rank(psi[i,], na.last="keep"))
  length_samples <- sum(!is.na(psi_rank))
  psi_low_tert <- which(psi_rank <= length_samples*1/3)
  psi_mid_tert <- which(psi_rank <= length_samples*2/3 & psi_rank > length_samples*1/3)
  psi_up_tert <- which(psi_rank > length_samples*2/3)
  
  num_supporting_samples <- c(num_supporting_samples, length_samples)
  
  psi[i, psi_up_tert] <- 1
  psi[i, psi_mid_tert] <- NA
  psi[i, psi_low_tert] <- 0
}

save(psi, gene_id, num_supporting_samples, file=paste0("../Input/", rbp, ".A5SS.PSI.discrete.Rdata"))

# A3SS event ####
load(file=paste0(rbp, ".A3SS-reads.Rdata"))

# Filter samples by requiring at least 10 reads of inclusion or exclusion 
# Filter that an event should be supported by at least 200 samples
# Filter splice targets with IQR(psi) > 0.1
incl_reads <- outReads[, seq(1,len,2)]
excl_reads <- outReads[, seq(2,len,2)]
total_reads <- incl_reads + excl_reads
denom <- incl_reads + 2*excl_reads
psi <- incl_reads / denom

excl_support <- rowSums(excl_reads >= 1) >= 1
table(excl_support)

total_support <- rowSums(total_reads >= 10) >= 1
table(total_support)

keep <- (excl_support + total_support) == 2
table(keep)

psi <- psi[keep,]
gene_id <- outName[keep,] # Make sure gene_id / outName stays consistent

iqr <- apply(psi, 1, function(x) IQR(x, na.rm=TRUE))
psi <- psi[iqr>0.10, ] # IQR of PSI > 0.10
gene_id <- gene_id[iqr>0.10, ]

num_supporting_samples <- numeric()

# save(psi, gene_id, file="A3SS.psi.Rdata")

for(i in 1:nrow(psi)){
  psi_rank <- as.numeric(rank(psi[i,], na.last="keep"))
  length_samples <- sum(!is.na(psi_rank))
  psi_low_tert <- which(psi_rank <= length_samples*1/3)
  psi_mid_tert <- which(psi_rank <= length_samples*2/3 & psi_rank > length_samples*1/3)
  psi_up_tert <- which(psi_rank > length_samples*2/3)
  
  num_supporting_samples <- c(num_supporting_samples, length_samples)
  
  psi[i, psi_up_tert] <- 1
  psi[i, psi_mid_tert] <- NA
  psi[i, psi_low_tert] <- 0
}

save(psi, gene_id, num_supporting_samples, file=paste0("../Input/", rbp, ".A3SS.PSI.discrete.Rdata"))

# RI event ####
load(file=paste0(rbp, ".RI-reads.Rdata"))

# Filter samples by requiring at least 10 reads of inclusion or exclusion 
# Filter that an event should be supported by at least 200 samples
# Filter splice targets with IQR(psi) > 0.1
incl_reads <- outReads[, seq(1,len,2)]
excl_reads <- outReads[, seq(2,len,2)]
total_reads <- incl_reads + excl_reads
denom <- incl_reads + 2*excl_reads
psi <- incl_reads / denom

excl_support <- rowSums(excl_reads >= 1) >= 1
table(excl_support)

total_support <- rowSums(total_reads >= 10) >= 1
table(total_support)

keep <- (excl_support + total_support) == 2
table(keep)

psi <- psi[keep,]
gene_id <- outName[keep,] # Make sure gene_id / outName stays consistent

iqr <- apply(psi, 1, function(x) IQR(x, na.rm=TRUE))
psi <- psi[iqr>0.10, ] # IQR of PSI > 0.10
gene_id <- gene_id[iqr>0.10, ]

num_supporting_samples <- numeric()

# save(psi, gene_id, file="RI.psi.Rdata")

for(i in 1:nrow(psi)){
  psi_rank <- as.numeric(rank(psi[i,], na.last="keep"))
  length_samples <- sum(!is.na(psi_rank))
  psi_low_tert <- which(psi_rank <= length_samples*1/3)
  psi_mid_tert <- which(psi_rank <= length_samples*2/3 & psi_rank > length_samples*1/3)
  psi_up_tert <- which(psi_rank > length_samples*2/3)
  
  num_supporting_samples <- c(num_supporting_samples, length_samples)
  
  psi[i, psi_up_tert] <- 1
  psi[i, psi_mid_tert] <- NA
  psi[i, psi_low_tert] <- 0
}

save(psi, gene_id, num_supporting_samples, file=paste0("../Input/", rbp, ".RI.PSI.discrete.Rdata"))

# MXE event ####
load(file=paste0(rbp, ".MXE-reads.Rdata"))

# Filter samples by requiring at least 10 reads of inclusion or exclusion 
# Filter that an event should be supported by at least 200 samples
# Filter splice targets with IQR(psi) > 0.1
incl_reads <- outReads[, seq(1,len,2)]
excl_reads <- outReads[, seq(2,len,2)]
total_reads <- incl_reads + excl_reads
denom <- incl_reads + 2*excl_reads
psi <- incl_reads / denom

excl_support <- rowSums(excl_reads >= 1) >= 1
table(excl_support)

total_support <- rowSums(total_reads >= 10) >= 1
table(total_support)

keep <- (excl_support + total_support) == 2
table(keep)

psi <- psi[keep,]
gene_id <- outName[keep,] # Make sure gene_id / outName stays consistent

iqr <- apply(psi, 1, function(x) IQR(x, na.rm=TRUE))
psi <- psi[iqr>0.10, ] # IQR of PSI > 0.10
gene_id <- gene_id[iqr>0.10, ]

num_supporting_samples <- numeric()

# save(psi, gene_id, file="MXE.psi.Rdata")

for(i in 1:nrow(psi)){
  psi_rank <- as.numeric(rank(psi[i,], na.last="keep"))
  length_samples <- sum(!is.na(psi_rank))
  psi_low_tert <- which(psi_rank <= length_samples*1/3)
  psi_mid_tert <- which(psi_rank <= length_samples*2/3 & psi_rank > length_samples*1/3)
  psi_up_tert <- which(psi_rank > length_samples*2/3)
  
  num_supporting_samples <- c(num_supporting_samples, length_samples)
  
  psi[i, psi_up_tert] <- 1
  psi[i, psi_mid_tert] <- NA
  psi[i, psi_low_tert] <- 0
}

save(psi, gene_id, num_supporting_samples, file=paste0("../Input/", rbp, ".MXE.PSI.discrete.Rdata"))

