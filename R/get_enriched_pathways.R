#' @export
get_enriched_pathways<-function(diff_genes,pathwayDB,expression_matrix_genes,pathwayDB_nodes,delta=0.05)
{
	enrichment_pvalue = c()
	enrichment_pvalue = sort(unlist(lapply(pathwayDB_nodes,function(x) get_enrichment_pvalue(x,diff_genes,expression_matrix_genes))))
	bins = list()
	if(length(which(p.adjust(enrichment_pvalue,method="BH") < delta)) > 0)
	{
		enrichment_pvalue = p.adjust(enrichment_pvalue,method="BH")
	}
	pathways_to_include = names(enrichment_pvalue)[enrichment_pvalue < delta]
	if(length(pathways_to_include) == 0){ return(NULL) }
	for(i in 1:length(pathways_to_include))
	{
		pathway_nodes = pathwayDB_nodes[[pathways_to_include[i]]]
		DEGs_in_pathway = intersect(names(diff_genes),pathway_nodes)
		bins[[i]] = DEGs_in_pathway
	}
	names(bins) = pathways_to_include
	return(bins)
}