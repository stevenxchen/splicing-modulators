# Splicing Modulators

We present a novel approach to interrogate large genomic datasets for modulators of the RBP-splicing relationship. 

In our model, we hypothesize that the splicing activity of an RBP can change with respect to putative modulators. RBP activity is estimated as the splicing levels of its target events, which are measured as percent spliced-in (PSI). Intuitively, if an RBP regulates the splicing outcome of a target event, we expect to observe a positive or negative correlation between the expression levels of the RBP and the PSI levels of the target events across multiple samples. A modulator candidate is selected if the correlation between the expression levels of a RBP and the PSI levels of a target event is dependent on the expression levels of a modulator. 

The strategy to evaluate modulator candidates integrates four categories of inputs: 
1) a large-scale RNAseq dataset 
2) an RBP of interest
3) corresponding RBP-binding sites
4) alternative splicing annotations. 

For an optional Mendelian randomization-based approach, matched genotype data are required.

The expression levels of the RBP and the PSI values of the target events are calculated from RNAseq data directly. The RBP-binding sites are derived from CLIP-seq data in the public domain. Splicing annotations are downloaded from GENCODE. The expression levels of the modulator candidates are imputed from genotype of the flanking regions.   

Please cite: