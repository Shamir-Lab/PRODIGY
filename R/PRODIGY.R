#' PRODIGY
#'
#' This function runs the PRODIGY algorithm for a single patient.
#' @param mutated_genes A vector of mutated genes to examine using PRODIGY.
#' @param expression_matrix A read count matrix with genes in rows and patients on columns. All genes must be contained in the global PPI network.
#' @param network The global PPI network. Columns describe the source protein, destination protein and interaction score respectively. The network is considered as undirected.
#' @param sample The sample label as appears in the SNV and expression matrices.
#' @param diff_genes A vector of the sample's differentially expressed genes (with gene names). All genes must be contained in the global PPI network.
#' @param alpha the penalty exponent.
#' @param pathway_list A list where each object is a 3 column data.table (src,dest,weight). Names correspond to pathway names
#' @param num_of_cores The number of CPU cores to be used by the influence scores calculation step.
#' @param sample_origins A vector that contains two optional values ("tumor","normal") corresponds to the tissues from which each column in expression_matrix was derived. This vector is utilized for differential expression analysis. If no vector is specified, the sample names of expression_matrix are assumed to be in TCGA format where last two digits correspond to sample type: "01"= solid tumor and "11"= normal.
#' @param write_results Should the results be written to text files?
#' @param results_folder Location for resulting influence matrices storage (if write_results = T)
#' @param beta Minimal fold-change threshold for declering gene as differentially expressed by DESeq (default = 0.2)
#' @param gamma FDR threshold for declering gene as differentially expressed by DESeq (default = 0.05)
#' @param delta FDR threshold for declering a pathway as statistically enriched for differentially expressed genes (default = 0.05)
#' @return A matrix of influence scores for every mutation and every enriched pathway.
#' @examples
#' # Load SNV+expression data from TCGA
#' data(COAD_SNV)
#' data(COAD_Expression)
#' # Load STRING network data 
#' data(STRING_network)
#' network = STRING_network
#' sample = intersect(colnames(expression_matrix),colnames(snv_matrix))[1]
#' # Identify sample origins (tumor or normal)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
#' res = PRODIGY<-function(snv_matrix,expression_matrix,network=network,sample,diff_genes=NULL,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins = sample_origins,beta=2,gamma=0.05,delta=0.05)
#' @references
#' Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 1-21 (2014).
#' Gabriele Sales, Enrica Calura and Chiara Romualdi, graphite: GRAPH Interaction from pathway Topological Environment (2017).
#' Gillespie, M., Vastrik, I., Eustachio, P. D., Schmidt, E. & Bono, B. De. Reactome: a knowledgebase of biological pathways. Nucleic Acids Res. 33, 428-432 (2005).
#' Schaefer, C. F. et al. PID: The pathway interaction database. Nucleic Acids Res. 37, 674-679 (2009).
#' Ogata, H. et al. KEGG: Kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 27, 29-34 (1999).
#' @export
PRODIGY<-function(mutated_genes,expression_matrix,network=NULL,sample,diff_genes=NULL,alpha=0.05,pathway_list = NULL,
			num_of_cores=1,sample_origins = NULL, write_results = F, results_folder = "./",
			beta=2,gamma=0.05,delta=0.05)
{
	#load needed R external packages
	libraries = c("DESeq2","igraph","ff","plyr","biomaRt","parallel","PCSF","graphite")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	#if no network is specified, the network derived from STRING is used as in the original publication
	if(is.null(network))
	{
		data(STRING_network)
		network = STRING_network
	}
	network[,"score"] = sapply(as.numeric(network[,"score"]),function(x) min(x,0.8))
	network[,"score"] = 1-as.numeric(network[,"score"])
	mutated_genes = mutated_genes[mutated_genes %in% unique(c(network[,1],network[,2]))]
	if(length(mutated_genes) < 1) {
	 	print("No mutated gene in large PPI network, aborting")
		return() 
	}
	original_network = network
	network = graph_from_data_frame(network,directed=F)
	if(is.null(pathway_list))
    {
		print("no pathway list. using Reactome as default")
		pathway_list = get_pathway_list_from_graphite(source = "reactome",minimal_number_of_nodes = 10,num_of_cores = num_of_cores)
	}
	#get differentially expressed genes list
	if(is.null(diff_genes))
	{
		expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
		#check if sample has expression data
		if(!(sample %in% colnames(expression_matrix)))
		{
			print("sample is missing expression information, aborting")
			return()
		}
		#if sample_origins = NULL we assume the matrices follow the TCGA convention (tumors suffix is "-01" and normals are "-11")
		if(is.null(sample_origins))
		{
			sample_origins = rep("tumor",ncol(expression_matrix))
			sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
		}
		if(length(which(sample_origins == "normal")) < 1)
		{	
			print("no normal samples, cannot perform differential expression analysis. aborting")
			return()
		}
		print("analyzing DEGs")	
		diff_genes = get_diff_expressed_genes(expression_matrix,sample,sample_origins,beta,gamma)
	}
    pathway_list_genes = setNames(lapply(pathway_list,function(x) get_pathway_genes_from_network(x)),names(pathway_list))
	#get enriched pathways. Bins contain the DEGs belonging to each pathway
	bins = get_enriched_pathways(diff_genes,pathwayDB = NULL,rownames(expression_matrix),pathwayDB_nodes = pathway_list_genes,delta)
	if(length(bins)==0){
		print("no enriched pathways, aborting")
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
		pathway_network = get_pathway_network(pathway_list[[pathway_name]],original_network)
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
	if(write_results)
	{
		write.table(Influence_matrix[-1,],file=paste(results_folder,sample,"_influence_scores.txt",sep="")
					,quote=F,col.names=T,row.names=T,sep="\t")
	}
	return(Influence_matrix[-1,])
}

#'' Return vector of gene names participate in the pathway
#'' @param pathway A 2 column data.table representing interactions
#'' return A vector of gene names 
#' @export
get_pathway_genes_from_network<-function(pathway)
{
    return(unique(c(pathway[,1],pathway[,2])))
}

#' Return a network induced by a given gene set
#' @param network A data.table with at least 2 columns. It is assumed that the first two columns contains the genes participating in the interaction
#' @param gene_set A vector of gene names.
#' return A data.table representing the network
#' @export
get_network_from_gene_set<-function(network,gene_set)
{
    return(network[network[,1] %in% gene_set & network[,2] %in% gene_set,])
}

#' Get list of data.table representing pathways using the graphite package. This could take a while so better to use more than a single core
#' @param source A string representing the graphite source to use. Default is Reactome
#' @param minimal_number_of_nodes An int, lower bound on the number of pathway genes
#' return A list of data.tables, each representing a single network
#' @export
get_pathway_list_from_graphite<-function(source = "reactome",minimal_number_of_nodes = 10,num_of_cores = 1)
{
	library(graphite)
	list_of_pathways = graphite::pathways("hsapiens",source)
	num_of_nodes = lapply(list_of_pathways,function(x) length(nodes(x)))
	pathway_names = names(list_of_pathways)[unlist(num_of_nodes) > minimal_number_of_nodes]
	if(num_of_cores == 1)
	{
		pathway_list = setNames(lapply(pathway_names,function(pathway_name) get_single_graphite_pathway(pathway_name,list_of_pathways)),pathway_names)	
	} else {
		pathway_list = setNames(mclapply(pathway_names,function(pathway_name) get_single_graphite_pathway(pathway_name,list_of_pathways),mc.cores=num_of_cores),pathway_names)	
	}
	pathway_list = pathway_list[sapply(pathway_list,length) > 0]
	return(pathway_list)
}

#' @export
get_single_graphite_pathway<-function(pathway_name,list_of_pathways)
{
	pathway_edges = tryCatch({
		return(graphite::edges(convertIdentifiers(list_of_pathways[[which(names(list_of_pathways) == pathway_name)]],"symbol"))[,c("src","dest")])
	},
	error = function(x){
		return(data.frame())
	},
	warning = function(x){
		return(data.frame())
	})
}
