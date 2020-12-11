# KTU
KTU: K-mer-based taxonomic clustering algorithm improves biological relevance in microbiome associated study

The 16S rDNA amplicon sequencing is widely implemented for microbiome associated studies. Microbiota feature picking algorithms for taxonomic identification and quantification are greatly renewing in recent years. The amplicon sequence variant (ASV) denoising algorithm of unbiased sequence picking replaces the OTU clustering methods. The ASV features can detect and distinguish the biological variations under the species OTU level (≧97% similarity). However, the quantification of single ASV among sequencing samples are sparse and less prevalent in the same biological groups. Here, we introduce a k-mer based, sequence alignment-free algorithm – “KTU” (K-mer Taxonomic Unit) for re-clustering ASVs into taxonomic units with more biological relevance.

The “KTU” algorithm was designed with four parts (k-mer frequency counting, k-mer frequency similarity measurement, k-mer feature partitioning, and generating KTU table) and conducted in the R environment. The k-mer frequency counting was conducted by tetranucleotide frequency of amplicons; 256 tetranucleotide compositions were then converted to a 0-to-1 proportion. The similarities of k-mer frequency among amplicons were measured by cosine similarity. The similarity matrix then was converted to the distance matrix for the subsequent step. The KTUs were detected from the cosine distance matrix by using partition around medoids (PAM) clustering algorithm; the iterative PAM-KTU detecting process found the optimal cluster numbers of KTUs. The final step of the KTU algorithm was aggregating ASVs into the KTUs and generating the KTU table.

We tested KTU algorithm with an open-access dataset from NCBI SRA. The SRA dataset was a wild flying squirrel fecal microbiota study that was associated with forest ambient temperature along with seasonal changes (Liu et al. 2018 from Microbial Ecology). We compared the microbiota composition associated to temperature changes using ASV dataset (r=0.52, p=0.022), KTU dataset (r=0.75, p<0.001), and the original results of paper (OTU clustering method; r=0.72, p< 0.001).

The 16S amplicon denoising method unbiased picks high-quality sequence variants but generates a sparse feature table for subsequent multivariate statistical analysis. We introduce a k-mer based re-clustering method for improving the environmental and biological relevancies of ASV based microbiome associated study effectively.


**Latest updates:** `r format(Sys.Date(), "%B %d, %Y")`

### Required packages:
`Biostrings` for parsing DNA sequence (fasta format) files
`coop` for measuring cosine similarity of k-mer frequency
`cluster` for conducting PAM/K-medoids clustering
`parallel` for parallel computing
`foreach` for parallel for-loop computing
`doParallel` for parallel computing
`digest` for generating MD5 hash


## Practice  
install package from github:
```
library(devtools)
install_github("poyuliu/KTU")
```

sourcing KTU R library  
~~source("http://www.lifescipy.net/RcodeDB/FN_KTU.R")~~  
```
library(KTU)
```

### Step 1 make k-mer DB (optional)  
Trim full length 16S DB to specific hypervariable regions  
Example: V3-V4 (341-805)  
Forward: **CCTACGGGNGGCWGCAG**  
Reverse: **GACTACHVGGGTATCTAATCC**  
output file: *SILVA_NR132_trimmed-sequence.fasta*  
```
trim.primer("/home/share/SILVA_NR132/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna",
            forseq="CCTACGGGNGGCWGCAG",
            revseq="GACTACHVGGGTATCTAATCC",
            output.name="SILVA_NR132")
```  
make and output DB
input: *trimmed fasta file*  
output: *SILVA_132_V3V4_DB.RDS & SILVA_132_V3V4_TX.RDS*  
```  
silva.db <- makektudb(input.fasta = "SILVA_NR132_trimmed-sequence.fasta",
                   input.taxa = "/home/share/SILVA_NR132/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt",
                   output.file = "SILVA_132_V3V4")

```

### Step 2 run Klustering (main algorithm)  
import ASV table (output from QIIME2 pipeline + DADA2 denoising + biom conversion) and remove original taxonomy column  
```
data <- read.delim("/home/poyuliu/NAS2/16S_SMI/16S_sarcopenia/feature-table_w_tax.txt")
data <- data[,-ncol(data)]
```  
*input representive-sequence fasta file*  
Parameters:  
  iterative PAM clustering steps: 5 and 10 (default:3 and 10)  
  CPU cores: 10 (default: 1)  
```
kluster <- klustering(repseq = "dna-sequences.fasta",
                          feature.table = data,
                          write.fasta = TRUE,step = c(5,10),cores = 10)
```  
Save results  
```
saveRDS(kluster,file="kluster.RDS")
```

### Step 3 assign k-mer taxonomy (optional)  
**similar to RDP-classifier/Naîve Bayesian algorithm**  
Read k-mer clustering result
```
kluster <- readRDS("kluster.RDS")
```  
run k-mer taxonomy assignment  
input DB: {name}_DB.RDS and {name}_TX.RDS  
input k-mer frequency table from previous result: *kluster$kmer.table*  

```
k.taxonomy <- kaxonomy(dbRDS = "SILVA_132_V3V4_DB.RDS",
                      taxaRDS = "SILVA_132_V3V4_TX.RDS",
                      kmer.table = kluster$kmer.table,
                      cores = 3)
```  
Save results  
```
write.table(k.taxonomy,file = "Kaxonomy.tsv",sep="\t")
```


## Example
[_Liu P.-Y., et al._ **Variations in Gut Microbiota of Siberian Flying Squirrels Correspond to Seasonal Phenological Changes in Their Hokkaido Subarctic Forest Ecosystem.** Microbial ecology 78, 223–231 (2019)](https://link.springer.com/article/10.1007/s00248-018-1278-x)


<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs00248-018-1278-x/MediaObjects/248_2018_1278_Fig2_HTML.png?as=webp" width="45%" /> <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs00248-018-1278-x/MediaObjects/248_2018_1278_Fig3_HTML.png?as=webp" width="45%" />

### Untreated ASV analyses
<img src="http://www.lifescipy.net/KTU/github/github_pic_01.png" width="60%" />
<img src="http://www.lifescipy.net/KTU/github/github_pic_02.png" width="60%" />
<img src="http://www.lifescipy.net/KTU/github/github_pic_03.png" width="60%" />

### KTU analyses
**Phyla composition**
<img src="http://www.lifescipy.net/KTU/github/github_pic_04.png" width="60%" />


**Beta diversity PCoA**
<img src="http://www.lifescipy.net/KTU/github/github_pic_05.png" width="60%" />


**Environmental factor correlation**
<img src="http://www.lifescipy.net/KTU/github/github_pic_06.png" width="60%" />
<img src="http://www.lifescipy.net/KTU/github/github_pic_07.png" width="60%" />


**Correlation heatmap**
<img src="http://www.lifescipy.net/KTU/github/github_pic_08.png" width="60%" />
