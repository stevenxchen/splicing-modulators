### Format and Analyze Output Data ####

rm(list=ls())
setwd("/Data/")
library(biomaRt)

setwd("Analyze")

splice.type = "SE" # Define splice type 

# load output data
load(paste0("/Data/SRSF1.Output.", splice.type, ".Rdata"))

imp_output <- output[complete.cases(output),]


imp_filter <- subset(imp_output, imp_output$gamma.pval < 0.05)

# Fisher's exact test
n.events <- length(unique(imp_filter$event))
x <- as.integer(0.05 * n.events)
f.table <-  matrix(c(x, as.integer(.05*n.events), n.events-x, as.integer(n.events-.05*n.events)), ncol=2)
fisher.test(f.table, alternative = "greater")
y <- numeric()
y <- c(y,fisher.test(f.table, alternative = "greater")$p.value)

top_imp_regs <- data.frame(Regulator_Ensembl = character(), Regulator = character(), Type = character(), Targets = integer(), stringsAsFactors = FALSE)
table_imp <- table(imp_filter$regulator)

for(j in 1:length(table_imp)){
  top_imp_regs[j,1] <- names(sort(table_imp, decreasing = T))[j]
  top_imp_regs[j,2] <- imp_filter$regulatorSymbol[match(top_imp_regs[j,1], imp_filter$regulator)]
  top_imp_regs[j,3] <- imp_filter$regulatorGeneType[match(top_imp_regs[j,1], imp_filter$regulator)]
  top_imp_regs[j,4] <- sort(table_imp, decreasing = T)[j]
}

top_imp_regs
write.csv(top_imp_regs, file = "temp.output.csv", row.names = FALSE)


### Plot histogram ####
library(ggplot2)
library(hrbrthemes)

ggplot(top_imp_regs, aes(x=Targets)) + 
  geom_histogram(binwidth=3, fill="#69b3a2", alpha=0.9) +
  # ggtitle("Histogram of SRSF1 Modulators and IR Targets") +
  theme_ipsum() +
  theme(
    #plot.title = element_text(size=15, hjust = "0.5")
    axis.text.y=element_text(size=7),
    axis.text.x=element_text(size=7)
  ) +
  ylab("Modulator Count") + xlab("IR Targets") +
  geom_vline(xintercept = 31, linetype="dashed", cex= 0.2)
 