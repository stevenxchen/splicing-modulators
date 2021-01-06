# Model to identify Modulators of RBP_Splicing ####
rm(list=ls())
setwd("/Data")

library(dplyr)

# Load Discrete Modulator Counts
load(file = "Modulator.discrete.Rdata")

# Load SRSF1 RBP Discrete Counts
load(file = "SRSF1.discrete.Rdata")
rbp_discrete <- as.integer(SRSF1.discrete)

# Define error function
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

# Define splice type
splice.type <- "SE"

### Model ####
load(file = paste0("SRSF1.", splice.type, ".PSI.discrete.Rdata"))

# Run code on all modulators and splice events ####
output <- data.frame(event=character(), GeneID=character(), geneSymbol=character(), 
                     RBP=character(), regulator=character(), samples=numeric(), gamma=numeric(), 
                     gamma.pval=numeric(), alpha=numeric(), alpha.pval=numeric(), 
                     beta=numeric(), beta.pval=numeric(), alpha.div.beta=numeric(), 
                     alphac=numeric(), alphaf=numeric(), betaf=numeric(),
                     stringsAsFactors = FALSE)
out_row <- 1

for(i in 1:nrow(Modulator.discrete)){ 
  reg_test <- as.numeric(Modulator.discrete[i,])
  for(j in 1:nrow(psi)){
    
    # Set the splice event
    psi_test <- as.numeric(psi[j,])
    
    # Temporary dataframe
    data_temp <- data.frame(mod=reg_test, rbp=rbp_discrete, psi=psi_test)
    
    pdat11 <- filter(data_temp, mod==1, rbp==1, !is.na(psi))
    p11 <- mean(pdat11$psi==1)
    
    pdat10 <- filter(data_temp, mod==1, rbp==0, !is.na(psi))
    p10 <- mean(pdat10$psi==1)
    
    pdat00 <- filter(data_temp, mod==0, rbp==0, !is.na(psi))
    p00 <- mean(pdat00$psi==1)
    
    pdat01 <- filter(data_temp, mod==0, rbp==1, !is.na(psi))
    p01 <- mean(pdat01$psi==1)
    
    gamma <- p11-p01-p10+p00
    betam <- p11-p01
    alpham <- p10-p00
    ab <- alpham/betam
    
    betaf <- p11-p10
    alphaf <- p01-p00
    alphac <- p00
    
    # p.value for alpham
    p.10.00 <- ( sum(pdat10$psi==1) + sum(pdat00$psi==1) ) / ( nrow(pdat10) + nrow(pdat00) )
    var.alpha <- p.10.00 * (1-p.10.00) * (1/nrow(pdat10) + 1/nrow(pdat00))
    pval.alpha <- erfc(abs(alpham/sqrt(2*var.alpha)))
    
    # p.value for betam
    p.11.01 <- ( sum(pdat11$psi==1) + sum(pdat01$psi==1) ) / ( nrow(pdat11) + nrow(pdat01) )
    var.beta <- p.11.01 * (1-p.11.01) * (1/nrow(pdat11) + 1/nrow(pdat01))
    pval.beta <- erfc(abs(betam/sqrt(2*var.beta)))
    
    # p.value for gamma
    p.gamma <- ( sum(pdat00$psi==1) + sum(pdat01$psi==1) + sum(pdat10$psi==1) + sum(pdat11$psi==1) ) /
      ( nrow(pdat00) + nrow(pdat01) + nrow(pdat10) + nrow(pdat11) )
    var.gamma <- p.gamma * (1-p.gamma) * (1/nrow(pdat00) + 1/nrow(pdat01) + 1/nrow(pdat10) + 1/nrow(pdat11))
    pval.gamma <- erfc(abs(gamma/sqrt(2*var.gamma)))
    
    # complete dataframe
    output[out_row,1] <- row.names(psi)[j]
    output[out_row,2] <- gene_id$name[j]
    output[out_row,3] <- gene_id$geneSymbol[j]
    output[out_row,4] <- "SRSF1"
    output[out_row,5] <- row.names(Modulator.discrete)[i]
    output[out_row,6] <- num_supporting_samples[j]
    output[out_row,7] <- gamma
    output[out_row,8] <- pval.gamma
    output[out_row,9] <- alpham
    output[out_row,10] <- pval.alpha
    output[out_row,11] <- betam
    output[out_row,12] <- pval.beta
    output[out_row,13] <- ab
    output[out_row,14] <- alphac
    output[out_row,15] <- alphaf
    output[out_row,16] <- betaf
    
    out_row <- out_row + 1
    
    print(c(i,j, splice.type))
  }
}

save(output, file = paste0("Output/SRSF1.Output.", splice.type, ".Rdata"))
