# KTU: K-mer-based taxonomic clustering algorithm improves biological relevance in microbiome associated study

The 16S rDNA amplicon sequencing is widely implemented for microbiome associated studies. Microbiota feature picking algorithms for taxonomic identification and quantification are greatly renewing in recent years. The amplicon sequence variant (ASV) denoising algorithm of unbiased sequence picking replaces the OTU clustering methods. The ASV features can detect and distinguish the biological variations under the species OTU level (≧97% similarity). However, the quantification of single ASV among sequencing samples are sparse and less prevalent in the same biological groups. Here, we introduce a k-mer based, sequence alignment-free algorithm – “KTU” (K-mer Taxonomic Unit) for re-clustering ASVs into taxonomic units with more biological relevance.

The “KTU” algorithm was designed with four parts (k-mer frequency counting, k-mer frequency similarity measurement, k-mer feature partitioning, and generating KTU table) and conducted in the R environment. The k-mer frequency counting was conducted by tetranucleotide frequency of amplicons; 256 tetranucleotide compositions were then converted to a 0-to-1 proportion. The similarities of k-mer frequency among amplicons were measured by cosine similarity. The similarity matrix then was converted to the distance matrix for the subsequent step. The KTUs were detected from the cosine distance matrix by using partition around medoids (PAM) clustering algorithm; the iterative PAM-KTU detecting process found the optimal cluster numbers of KTUs. The final step of the KTU algorithm was aggregating ASVs into the KTUs and generating the KTU table.

We tested KTU algorithm with an open-access dataset from NCBI SRA. The SRA dataset was a wild flying squirrel fecal microbiota study that was associated with forest ambient temperature along with seasonal changes (Liu et al. 2018 from Microbial Ecology). We compared the microbiota composition associated to temperature changes using ASV dataset (r=0.52, p=0.022), KTU dataset (r=0.75, p<0.001), and the original results of paper (OTU clustering method; r=0.72, p< 0.001).

The 16S amplicon denoising method unbiased picks high-quality sequence variants but generates a sparse feature table for subsequent multivariate statistical analysis. We introduce a k-mer based re-clustering method for improving the environmental and biological relevancies of ASV based microbiome associated study effectively.

##  **More updated information**
This package has been published in _Methods in Ecology and Evolution_. The updated tutorial (with format description) and example files are available, please find them from the following links:
https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13758  
https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13758&file=mee313758-sup-0001-Supinfo.pdf  

**How to cite the KTU algorithm and paper:**  
Liu, P.-Y., Yang, S.-H., & Yang, S.-Y. (2022). KTU: K-mer Taxonomic Units improve the biological relevance of amplicon sequence variant microbiota data. _Methods in Ecology and Evolution_, 13, 560– 568. https://doi.org/10.1111/2041-210X.13758

**Example dataset in the paper: microbiome of sourdough starters**  
This repository includes a feature table (tab-delimited text file, ASV hash ID in the first column) and representative sequences (fasta format): https://github.com/poyuliu/KTU-validation/tree/main/sourdough
_Note that the feature table should include the first column with ASV IDs._

**Please find out the tutorial on wiki**  
https://github.com/poyuliu/KTU/wiki

Many thanks for [@csmiguel](https://github.com/csmiguel) providing codes for a workflow of KTU integration with dada2.  
https://github.com/poyuliu/KTU/wiki/2.-KTU-integration-with-dada2-workflow
