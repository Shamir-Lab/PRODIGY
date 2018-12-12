#' PRODIGY
#'
#' This function runs the PRODIGY algorithm for a single patient.
#' @param snp_matrix A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
#' @param expression_matrix A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
#' @param network The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
#' @param sample The sample label as appears in the SNP and expression matrices.
#' @param diff_genes A vector of the sample's differentially expressed genes (with gene names). All genes must be contained in the global PPI network.
#' @param alpha the penalty exponent.
#' @param pathwayDB The pathway DB name from which curated pathways are taken. Could be one of three built in reservoirs ("reactome","kegg","nci").
#' @param num_of_cores The number of CPU cores to be used by the influence scores calculation step.
#' @param sample_origins A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
#' @return A matrix of influence scores for every mutation and every enriched pathway.
#' @examples
#' # Load SNP+expression data from TCGA
#' data(COAD_SNP)
#' data(COAD_Expression)
#' # Load STRING network data 
#' data(STRING_network)
#' network = STRING_network
#' sample = intersect(colnames(expression_matrix),colnames(snp_matrix))[1]
# Identify sample origins (tumor or normal)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
#' res = PRODIGY<-function(snp_matrix,expression_matrix,network=network,sample,diff_genes=NULL,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins = sample_origins)
#' @references
#' Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1-21 (2014).
#' Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).
#' Gillespie, M., Vastrik, I., Eustachio, P. D., Schmidt, E. & Bono, B. De. Reactome?: a knowledgebase of biological pathways. Nucleic Acids Res. 33, 428-432 (2005).
#' Schaefer, C. F. et al. PID: The pathway interaction database. Nucleic Acids Res. 37, 674-679 (2009).
#' Ogata, H. et al. KEGG: Kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 27, 29-34 (1999).
#' @export
PRODIGY<-function(snp_matrix,expression_matrix,network=NULL,sample,diff_genes=NULL,alpha=0.05,pathwayDB="reactome",
			num_of_cores=1,sample_origins = NULL)
{
	#load needed R external packages
	libraries = c("DESeq2","igraph","ff","plyr","biomaRt","parallel","PCSF")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	#check if sample has both SNP and expression data
	if(!(sample %in% colnames(snp_matrix) & sample %in% colnames(expression_matrix)))
	{
		print("sample is missing SNP or expression information, aborting")
		return()
	}
	#if no network is specified, the network derived from STRING is used as in the original publication
	if(is.null(network))
	{
		data(STRING_network.RData)
		network = STRING_network
	}
	network[,"score"] = min(as.numeric(network[,"score"]),0.8)
	network[,"score"] = 1-as.numeric(network[,"score"])
	expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
	snp_matrix = snp_matrix[which(rownames(snp_matrix) %in% unique(c(network[,1],network[,2]))),]
	original_network = network
	network = graph_from_data_frame(network,directed=F)
	#if sample_origins = NULL we assume the matrices follow the TCGA convention (tumors prefix is "-01" and normals are "-11")
	if(is.null(sample_origins))
	{
		sample_origins = rep("tumor",ncol(expression_matrix))
		sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
	}
	# if pathwayDB_nodes is not spacified, we use predefined lists from Reactome, KEGG or NCI PID
	# all pathways were acquired using the "graphite" R package (Gabriele Sales, Enrica Calura and Chiara Romualdi (2017))	
	data(pathwayDB_nodes)
	if(pathwayDB == "reactome"){	pathwayDB_nodes = pathwayDB_nodes[["reactome"]] 
	}else if(pathwayDB == "kegg") { pathwayDB_nodes = pathwayDB_nodes[["kegg"]]
	}else if(pathwayDB == "nci"){	pathwayDB_nodes = pathwayDB_nodes[["nci"]] 
	} else { 
		print("no pathwayDB selected or defined. aborting")
		return()
	}
	#get differentially expressed genes list
	if(is.null(diff_genes))
	{
		if(length(which(sample_origins == "normal")) < 1)
		{	
			print("no normal samples, cannot perform differential expression analysis. aborting")
			return()
		}
		print("analyzing DEGs")	
		diff_genes = get_diff_expressed_genes(expression_matrix,sample,sample_origins)
	}
	#get enriched pathways. Bins contain the DEGs belonging to each pathway
	bins = get_enriched_pathways(diff_genes,pathwayDB,rownames(expression_matrix),pathwayDB_nodes)
	if(length(bins)==0){
		print("no enriched pathways, aborting")
		return() 
	}
	mutated_genes = rownames(snp_matrix)[which(snp_matrix[,sample]==1)]
	if(length(mutated_genes) < 2) {
	 	print("no mutations to rank, aborting")
		return() 
	}
	print(paste("checking",length(mutated_genes),"mutations and",length(bins),"pathways"))
	Influence_matrix = matrix(ncol = length(mutated_genes))
	colnames(Influence_matrix) = mutated_genes
	Fraction_of_DEGs_matrix = matrix(ncol = length(mutated_genes))
	colnames(Fraction_of_DEGs_matrix) = mutated_genes
	#Main part of the algorithm: calculate influence scores for every mutation and every 
	#deregulated pathway
	for(i in 1:length(bins))
	{
		pathway_name = names(bins)[i]
		pathway_id = pathway_name
		seed = bins[[i]]
		pathway_network = get_pathway_network(pathwayDB,pathway_name,original_network)
		if(is.null(nrow(pathway_network))){ next }
		#prize nodes are DEGs belong to the pathway
		pathway_prizes = diff_genes[bins[[pathway_name]]]
		max_score = sum(diff_genes[seed])
		res=matrix(unlist(mclapply(mutated_genes,function(x) 
			run_single_PCSF(pathway_network,x,original_network,pathway_prizes,alpha),mc.cores = num_of_cores)),
			ncol=2,byrow=T)
		if(nrow(res) < length(mutated_genes)){ next }
		influence_scores = res[,1]/max_score
		fraction_of_DEGs_values = res[,2]/max_score
		Influence_matrix = rbind(Influence_matrix,influence_scores)
		rownames(Influence_matrix)[nrow(Influence_matrix )] = pathway_name 
		Fraction_of_DEGs_matrix = rbind(Fraction_of_DEGs_matrix ,fraction_of_DEGs_values)
		rownames(Fraction_of_DEGs_matrix)[nrow(Fraction_of_DEGs_matrix )] = pathway_name 
		gc()
	}
	return(Influence_matrix[-1,])
}
