<!-- toc -->

דצמבר 12, 2018

# DESCRIPTION

```
Package: PRODIGY
Title: Personalized prioritization of driver genes
Version: 1.0
Authors: Gal Dinstag, Ron Shamir
Reference: Dinstag, G. & Shamir R. PRODIGY: personalized prioritization of driver genes. bioRxiv (2018)  
Maintainer: Gal Dinstag <galdinstag@mail.tau.ac.il>
Description: PRODIGY analyzes the expression and SNV profiles of a patient along with data on known pathways and protein-protein interactions. PRODIGY quantifies the impact of each mutated gene on every transcriptionally deregulated pathway using the prize collecting Steiner tree model. Mutated genes are ranked by their aggregated impact on all deregulated pathways.
Depends: R (>= 3.4.3)
License: GPL
Encoding: UTF-8
LazyData: true
RoxygenNote: 6.1.1```


# `analyze_PRODIGY_results`

analyze_PRODIGY_results


## Description

This function analyzes PRODIGY's results (influence score matrices) to rank individual driver genes for a group of patients.


## Usage

```r
analyze_PRODIGY_results(all_patients_scores)
```


## Arguments

Argument      |Description
------------- |----------------
`all_patients_scores`     |     Either a list of influence scores matrices or a single matrix received by the PRODIGY algorithm.


## Value

A list of ranked genes per sample.


## Examples

```r
results = analyze_PRODIGY_results(all_patients_scores) #' @references
```


# `get_DEGs`

get_DEGs


## Description

This function calculates differentially expressed genes (DEGs) using DESeq2 for a group of patients.  DEGs calculation can take a while so it is recommended to make this analysis as a pre-process


## Usage

```r
get_DEGs(expression_matrix, samples, sample_origins = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`expression_matrix`     |     A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
`samples`     |     The sample labels as they appear in the expression matrix.
`sample_origins`     |     a vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.


## Value

A named list of DEGs per sample.


## References

Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1–21 (2014).


## Examples

```r
data(COAD_Expression)
sample_origins = rep("tumor",ncol(expression_matrix))
sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"
expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
```


# `PRODIGY`

PRODIGY


## Description

This function runs the PRODIGY algorithm for a single patient.


## Usage

```r
PRODIGY(snp_matrix, expression_matrix, network = NULL, sample,
  diff_genes = NULL, alpha = 0.05, pathwayDB = "reactome",
  num_of_cores = 1, sample_origins = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`snp_matrix`     |     A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
`expression_matrix`     |     A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
`network`     |     The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
`sample`     |     The sample label as appears in the SNP and expression matrices.
`diff_genes`     |     A vector of the sample's differentially expressed genes (with gene names). All genes must be contained in the global PPI network.
`alpha`     |     the penalty exponent.
`pathwayDB`     |     The pathway DB name from which curated pathways are taken. Could be one of three built in reservoirs ("reactome","kegg","nci").
`num_of_cores`     |     The number of CPU cores to be used by the influence scores calculation step.
`sample_origins`     |     A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.


## Value

A matrix of influence scores for every mutation and every enriched pathway.


## References

Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1–21 (2014).
 Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).
 Gillespie, M., Vastrik, I., Eustachio, P. D., Schmidt, E. & Bono, B. De. Reactome?: a knowledgebase of biological pathways. Nucleic Acids Res. 33, 428–432 (2005).
 Schaefer, C. F. et al. PID: The pathway interaction database. Nucleic Acids Res. 37, 674–679 (2009).
 Ogata, H. et al. KEGG: Kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 27, 29–34 (1999).


## Examples

```r
# Load SNP+expression data from TCGA
data(COAD_SNP)
data(COAD_Expression)
# Load STRING network data
data(STRING_network)
network = STRING_network
sample = intersect(colnames(expression_matrix),colnames(snp_matrix))[1]
sample_origins = rep("tumor",ncol(expression_matrix))
sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"
res = PRODIGY<-function(snp_matrix,expression_matrix,network=network,sample,diff_genes=NULL,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins = sample_origins)
```


# `PRODIGY_cohort`

PRODIGY_cohort


## Description

This is a wrapper that runs PRODIGY for a cohort of patients.


## Usage

```r
PRODIGY_cohort(snp_matrix, expression_matrix, network = NULL,
  samples = NULL, DEGs = NULL, results_folder = "./", alpha = 0.05,
  pathwayDB = "reactome", num_of_cores = 1, sample_origins = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`snp_matrix`     |     A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
`expression_matrix`     |     A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
`network`     |     The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
`DEGs`     |     Named list of differentially expressed genes for every sample. All genes must be contained in the global PPI network.
`alpha`     |     The penalty exponent.
`pathwayDB`     |     The pathway DB name from which curated pathways are taken. Could be one of three built in reservoirs ("reactome","kegg","nci").
`num_of_cores`     |     The number of CPU cores to be used by the influence scores calculation step.
`sample_origins`     |     A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
`sample`     |     The sample labels as appears in the SNP and expression matrices.


## Value

A list of influence scores matrices.


## References

Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1–21 (2014).
 Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).


## Examples

```r
data(COAD_SNP)
data(COAD_Expression)
# Load STRING network data
data(STRING_network)
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
all_patients_scores = PRODIGY_cohort(snp_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins)
```


