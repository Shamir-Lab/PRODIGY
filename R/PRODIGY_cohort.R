#' PRODIGY_cohort
#'
#' This is a wrapper that runs PRODIGY for a cohort of patients. It could take a long time to finish a group of patints, hence it is
#' advised to output PRODIGY's results into files using the write_results parameter. It is also advised to use as much cores as possible
#' using the num_of_cores parameter
#' @param snv_matrix A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
#' @param expression_matrix A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
#' @param network The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
#' @param sample The sample labels as appears in the SNV and expression matrices.
#' @param DEGs Named list of differentially expressed genes for every sample. All genes must be contained in the global PPI network.
#' @param alpha The penalty exponent.
#' @param pathway_list A list where each object is a 3 column data.table (src,dest,weight). Names correspond to pathway names
#' @param num_of_cores The number of CPU cores to be used by the influence scores calculation step.
#' @param sample_origins A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
#' @param write_results Should the results be written to text files?
#' @param results_folder Location for resulting influence matrices storage (if write_results = T) 
#' @param beta Minimal fold-change threshold for declering gene as differentially expressed by DESeq (default = 0.2)
#' @param gamma FDR threshold for declering gene as differentially expressed by DESeq (default = 0.05)
#' @param delta FDR threshold for declering a pathway as statistically enriched for differentially expressed genes (default = 0.05)
#' @param mutation_list An alternative to snv_matrix, one can simply provide a list of named vectors with genes to be analyzed by PRODIGY
#' @return A list of influence scores matrices.
#' @examples
#' data(COAD_SNV)
#' data(COAD_Expression)
#' # Load STRING network data 
#' data(STRING_network)
#' network = STRING_network
#' # Take samples for which SNV and expression is available 
#' samples = intersect(colnames(expression_matrix),colnames(snv_matrix))[1:5]
#' # Get differentially expressed genes (DEGs) for all samples
#' expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
#' DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
#' # Identify sample origins (tumor or normal)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
#' # Run PRODIGY
#' all_patients_scores = PRODIGY_cohort(snv_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins,
#' write_results = F, results_folder = "./",beta=2,gama=0.05,delta=0.05)
#' @references
#' Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1-21 (2014).
#' Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).
#' @export
PRODIGY_cohort<-function(snv_matrix = NULL,expression_matrix,network=NULL,samples=NULL,DEGs=NULL,alpha=0.05,pathway_list=NULL,
			num_of_cores=1,sample_origins = NULL, write_results = F, results_folder = "./",beta=2,gamma=0.05,delta=0.05,mutation_list = NULL)
{
	#load needed R external packages
	libraries = c("DESeq2","igraph","ff","plyr","biomaRt","parallel","PCSF","graphite")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	if(is.null(samples))
	{
		print("no samples, aborting")
		return()
	}
	all_patients_scores = list()
	if(is.null(sample_origins))
	{
		sample_origins = rep("tumor",ncol(expression_matrix))
		sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
	}
	if(is.null(network))
	{
		data(STRING_network)
		network = STRING_network
	}
	#run PRODIGY for all samples
	for(sample in samples)
	{
	   print(sample)
	   if(!is.null(DEGs) & (sample %in% names(DEGs))) 
	   {	 
			diff_genes = DEGs[[sample]] 
	   } else { diff_genes = NULL }
			if(!is.null(mutation_list))
			{
				sample_mutations = mutation_list[[sample]]
			} else {
				sample_mutations = names(snv_matrix[snv_matrix[,sample] == 1,sample])
			}
			all_patients_scores[[sample]] = PRODIGY(sample_mutations,expression_matrix,network,sample,diff_genes,alpha=alpha,pathway_list=pathway_list,
	                                            num_of_cores=num_of_cores,sample_origins = sample_origins,write_results = write_results,
												results_folder = results_folder,
	                                            beta = beta, gamma = gamma, delta = delta)
	}
	return(all_patients_scores)
}