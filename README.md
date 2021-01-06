# splicing-modulators

We present a novel approach to interrogate large genomic datasets for modulators of the RBP-splicing relationship. 

In our model, we hypothesize that the splicing activity of an RBP can change with respect to putative modulators. RBP activity is estimated as the splicing levels of its target events, which are measured as percent spliced-in (PSI). Intuitively, if an RBP regulates the splicing outcome of a target event, we expect to observe a positive or negative correlation between the expression levels of the RBP and the PSI levels of the target events across multiple samples. A modulator candidate is selected if the correlation between the expression levels of a RBP and the PSI levels of a target event is dependent on the expression levels of a modulator. 

The strategy to evaluate modulator candidates integrates four categories of inputs: a large-scale RNAseq dataset paired with genotype data, an RBP of interest, corresponding RBP-binding sites, and alternative splicing annotations. The expression levels of the RBP and the PSI values of the target events are calculated from RNAseq data directly. The RBP-binding sites are derived from CLIP-seq data in the public domain. Splicing annotations are downloaded from GENCODE. The expression levels of the modulator candidates are imputed from genotype of the flanking regions.   

The modulator candidate is identified using a generalized linear model: 

Yt=𝛽o+𝛽1Xr+𝛽2Xm+𝛽3XrXm+𝜀

where:
Xr is the gene expression of the RBP, 
Xm  is the imputed expression level of a modulator candidate, and 
Yt  is the PSI value of one of the RBP target event. 
Non-zero outcomes of 𝛽3 represent interactions between the imputed gene expression level of modulator and the gene expression level of RBP on the given splicing event.  

Paper is available at:

Please cite: