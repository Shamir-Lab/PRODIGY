#' @export
get_pathway_network<-function(pathwayDB,pathway_name,original_network)
{
	list_of_pathways = pathways("hsapiens", pathwayDB)
	pathway_graph = graphite::edges(convertIdentifiers(list_of_pathways[[which(names(list_of_pathways) == pathway_name)]],"symbol"))[,c(1,2)]
	#delete duplicated edges
	matches = match_df(data.frame(src=original_network[,1],dest=original_network[,2]),rbind(pathway_graph[,c(1,2)],pathway_graph[,c(2,1)]))
	if(nrow(matches) > 0) { 
		original_network = original_network[-as.numeric(rownames(matches)),] 
	}
	#add the edges to the background network
	#all edges have constant score of 0.1
	original_network = rbind(original_network,matrix(c(pathway_graph[,1],pathway_graph[,2],
					rep(0.1,nrow(pathway_graph))),ncol=3,nrow=nrow(pathway_graph)))
	#delete self-loops
	original_network = original_network[which(original_network[,1] != original_network[,2]),]
	#use only edgeg where one of the ends (or both) is in the pathway. 
	pathway_nodes = unique(c(as.character(pathway_graph[,1]),as.character(pathway_graph[,2])))
	pathway_network = original_network[which(original_network[,1] %in% pathway_nodes | original_network[,2] %in% pathway_nodes),]
	return(pathway_network)
}