#' @export
calculate_precision<-function(ranked_list,phenotype_genes)
{
	precision_vector = c()
	for(i in 1:length(ranked_list))
	{
		precision_vector = c(precision_vector,(length(intersect(ranked_list[1:i],phenotype_genes)))/
						i)
	}
	precision_vector[c(length(ranked_list):100)] = precision_vector[length(ranked_list)]
	return(precision_vector)
}

