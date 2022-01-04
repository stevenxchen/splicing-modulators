# Splicing Modulators

We present a novel approach to interrogate large genomic datasets for modulators of the RBP-splicing relationship.

In our model, we hypothesize that the splicing activity of an RBP can change with respect to putative modulators. RBP activity is estimated as the splicing levels of its target events, which are measured as percent spliced-in (PSI). Intuitively, if an RBP regulates the splicing outcome of a target event, we expect to observe a positive or negative correlation between the expression levels of the RBP and the PSI levels of the target events across multiple samples. A modulator candidate is selected if the correlation between the expression levels of a RBP and the PSI levels of a target event is dependent on the expression levels of a modulator.

The strategy to evaluate modulator candidates integrates four categories of inputs:
1) a large-scale RNAseq dataset
2) an RBP of interest
3) corresponding RBP-binding sites
4) alternative splicing annotations.

For an optional Mendelian randomization-based approach, a database of transcrption prediction weights appropriate to the sample type is required. A repository of weights can be found at https://predictdb.org/. Imputation of gene expression can be accomplished using PrediXcan, now part of the MetaXcan toolset https://github.com/hakyimlab/MetaXcan.

To run the software, you will need a RNAseq dataset with a large number of samples and a computing environment capable of processing the data. In this study, we used the ROSMAP and CommonMind Consortium datasets. Both are available on the AD Knowledge Portal: https://adknowledgeportal.synapse.org/.

To prepare the data, code in the "splicing_annotations" folder will count the number of reads for each splicing event. It requires a file input of splicing annotations (GeneID, chr, strand, start location, end location) for each splicing type (SE, MXE, A5SS, A3SS, RI). These inputs can be extracted from the comprehensive gene annotation GTF file for the corresponding genome reference from GENCODE https://www.gencodegenes.org/.

Next, RNA-binding protein sites from CLIP data are identified and overlapped with splicing annotations. Public domain CLIP data may be found on ENCODE https://www.encodeproject.org/.

Finally, all inputs are used to run the program in the "modulator_selection" folder. Data are initialized into discrete variables and then run accordingly using our algorithm.


For more information, please visit the scientific paper, published at

Please cite:
