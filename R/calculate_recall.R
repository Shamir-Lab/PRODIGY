#' @export
calculate_recall<-function(ranked_list,phenotype_genes)
{
	recall_vector = c()
	for(i in 1:length(ranked_list))
	{
		recall_vector = c(recall_vector,(length(intersect(ranked_list[1:i],phenotype_genes)))/
						length(phenotype_genes))
	}
	recall_vector[c(length(ranked_list):100)] = recall_vector[length(ranked_list)]
	return(recall_vector)
}