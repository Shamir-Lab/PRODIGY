
# PRODIGY
This R package prioritize driver genes for individual cancer patients.

The details of the method are described in
Dinstag G. & Shamir R. PRODIGY: personalized prioritization of driver genes. Bioinformatics (2019), 
https://academic.oup.com/bioinformatics/article/36/6/1831/5612092

## Package installation
```r
library(devtools)
install_github("Shamir-Lab/PRODIGY")
```

PRODIGY was developed using and is dependent on the following packages (minimal version required is specified):

- MASS_7.3-50
- DESeq2_1.16.1
- igraph_1.2.2
- graphite_1.22
- ff_2.2-14
- plyr_1.8.4
- biomaRt_2.32.1
- PCSF_0.99.1
- mixtools_1.1.0
- ggplot2_3.0.0
- cowplot_0.9.3

## Simple run example
```r
library(PRODIGY)
# Load SNP+expression data derived from TCGA
data(COAD_SNV)
data(COAD_Expression)
# Load STRING network data 
data(STRING_network)
network = STRING_network
# Take samples for which SNP and expression is available 
samples = intersect(colnames(expression_matrix),colnames(snv_matrix))[1:5]
# Get differentially expressed genes (DEGs) for all samples
expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL,beta=2,gamma=0.05)
# Identify sample origins (tumor or normal)
sample_origins = rep("tumor",ncol(expression_matrix))
sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
list_of_pathways = get_pathway_list_from_graphite(source = "reactome",minimal_number_of_nodes = 10,num_of_cores = 20)
# Run PRODIGY
all_patients_scores = PRODIGY_cohort(snv_matrix = snv_matrix,expression_matrix = expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,
 			pathway_list=list_of_pathways,num_of_cores=30,sample_origins=sample_origins,write_results = F, results_folder = "./",
 			beta=2,gamma=0.05,delta=0.05)
# Get driver gene rankings for all samples 
results = analyze_PRODIGY_results(all_patients_scores) 
```
