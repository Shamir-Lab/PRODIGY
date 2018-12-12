# PRODIGY
This R package prioritize driver genes for individual cancer patients.

The details of the method are described in
Dinstag, G. & Shamir R. PRODIGY: personalized prioritization of driver genes. bioRxiv (2018)
https://www.biorxiv.org/content/early/2018/10/30/456723

## Package installation
```r
library(devtools)
install_github("Shamir-Lab/PRODIGY")
```

Needed packages and minimal versions:

- MASS_7.3-50
- DESeq2_1.16.1
- igraph_1.2.2
- ff_2.2-14
- plyr_1.8.4
- biomaRt_2.32.1
- PCSF_0.99.1
- mixtools_1.1.0

## Simple run example
```r
library(PRODIGY)
# Load SNP+expression data derived from TCGA
load("data/COAD_SNP.RData")
load("data/COAD_Expression.RData")
# Load STRING network data 
load("data/STRING_network.RData")
network = STRING_network
# Take samples for which SNP and expression is available 
samples = intersect(colnames(expression_matrix),colnames(snp_matrix))[1:5]
# Get differentially expressed genes (DEGs) for all samples
expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
# Identify sample origins (tumor or normal)
sample_origins = rep("tumor",ncol(expression_matrix))
sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
# Run PRODIGY
all_patients_scores = PRODIGY_cohort(snp_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,
			pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins)
# Get driver gene rankings for all samples 
results = analyze_PRODIGY_results(all_patients_scores) 
```