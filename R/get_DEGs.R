#' get_DEGs
#'
#' This function calculates differentially expressed genes (DEGs) using DESeq2 for a group of patients.  DEGs calculation can take a while so it is recommended to make this analysis as a pre-process
#' @param expression_matrix A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
#' @param samples The sample labels as they appear in the expression matrix.
#' @param sample_origins a vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
#' @param beta Minimal fold-change threshold for declering gene as differentially expressed by DESeq (default = 0.2)
#' @param gamma FDR threshold for declering gene as differentially expressed by DESeq (default = 0.05)
#' @return A named list of DEGs per sample.
#' @examples
#' data(COAD_Expression)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"
#' expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
#' DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
#' @references
#' Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1-21 (2014).
#' @export
get_DEGs<-function(expression_matrix,samples,sample_origins=NULL,beta=2,gamma=0.05)
{
	libraries = c("DESeq2")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	if(is.null(sample_origins))
	{
		sample_origins = rep("tumor",ncol(expression_matrix))
		sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
		if(length(which(sample_origins == "normal")) < 1)
		{	
			print("no normal samples, cannot perform differential expression analysis. aborting")
			return()
		}
	}
	DEGs = list()
	for(sample in samples)
	{
		diff_genes = get_diff_expressed_genes(expression_matrix,sample,sample_origins,beta,gamma)
		DEGs[[sample]] = diff_genes 
	}
	return(DEGs)
}