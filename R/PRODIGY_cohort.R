#' PRODIGY_cohort
#'
#' This is a wrapper that runs PRODIGY for a cohort of patients.
#' @param snp_matrix A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
#' @param expression_matrix A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
#' @param network The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
#' @param sample The sample labels as appears in the SNP and expression matrices.
#' @param DEGs Named list of differentially expressed genes for every sample. All genes must be contained in the global PPI network.
#' @param alpha The penalty exponent.
#' @param pathwayDB The pathway DB name from which curated pathways are taken. Could be one of three built in reservoirs ("reactome","kegg","nci").
#' @param num_of_cores The number of CPU cores to be used by the influence scores calculation step.
#' @param sample_origins A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
#' @return A list of influence scores matrices.
#' @examples
#' load("data/COAD_SNP.RData")
#' load("data/COAD_Expression.RData")
#' # Load STRING network data 
#' load("data/STRING_network.RData")
#' network = STRING_network
#' # Take samples for which SNP and expression is available 
#' samples = intersect(colnames(expression_matrix),colnames(snp_matrix))[1:5]
#' # Get differentially expressed genes (DEGs) for all samples
#' expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
#' DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
#' # Identify sample origins (tumor or normal)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
#' # Run PRODIGY
#' all_patients_scores = PRODIGY_cohort(snp_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins)
#' @references
#' Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1–21 (2014).
#' Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).

PRODIGY_cohort<-function(snp_matrix,expression_matrix,network=NULL,samples=NULL,DEGs=NULL,results_folder = "./",alpha=0.05,pathwayDB="reactome",
			num_of_cores=1,sample_origins = NULL)
{
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
	# we use predefined lists of pathways from Reactome, KEGG or NCI PID
	# all pathways were acquired using the "graphite" R package (Gabriele Sales, Enrica Calura and Chiara Romualdi (2017))
	if(pathwayDB == "reactome"){	load("data/Reactome.RData") 
	}else if(pathwayDB == "kegg") { load("data/Kegg.RData")
	}else if(pathwayDB == "nci"){	load("data/NCI_PID.RData") 
	} else { 
		print("no pathwayDB selected. aborting")
		return()
	}
	if(is.null(network))
	{
		load("data/STRING_network.RData")
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
		all_patients_scores[[sample]] = PRODIGY(snp_matrix,expression_matrix,network,sample,diff_genes,alpha=alpha,pathwayDB=pathwayDB,
			num_of_cores=num_of_cores,sample_origins = sample_origins)
	}
	return(all_patients_scores)
}