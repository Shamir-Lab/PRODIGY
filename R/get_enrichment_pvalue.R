#' @export
get_enrichment_pvalue<-function(pathway_nodes,diff_genes,expression_matrix_genes)
{
	target_set = intersect(pathway_nodes,expression_matrix_genes)
	if(length(target_set) == 0)
	{
		return(1) 
	}
	background_set = setdiff(expression_matrix_genes,target_set)
	drawn = names(diff_genes)
	special_drawn = intersect(target_set,drawn)
	if(length(special_drawn) == 0)
	{
		return(1) 
	}
	return(phyper(length(special_drawn),length(target_set),length(background_set),length(drawn), lower.tail = F))	
}