#' @export
run_single_PCSF<-function(pathway_network,driver_gene,original_network,pathway_prizes,alpha)
{
	#add edges that touch the driver and add negative prizes for these nodes
	pathway_nodes = unique(c(pathway_network[,1],pathway_network[,2]))
	lines_to_add = which(original_network[,1]==driver_gene | original_network[,2]==driver_gene)
	#dont add edges that already touch pathway nodes because we already added them
	pathway_network = rbind(pathway_network,matrix(c(original_network[lines_to_add,1],
				original_network[lines_to_add,2],as.numeric(original_network[lines_to_add,3])),nrow=length(lines_to_add),ncol=3))
	if(length(which(duplicated(pathway_network))) > 0) { 
		pathway_network = pathway_network[-which(duplicated(pathway_network)),] 
	}
	pathway_network_igraph = graph_from_data_frame(as.data.frame(pathway_network[,c(1,2)]),directed=F)
	#if the putative driver cannot reach any gene in the pathway there is no need to calculate influence score
	if(all(is.infinite(distances(pathway_network_igraph,v = driver_gene,to = pathway_nodes))))
	{
		return(c(0,0))
	}
	E(pathway_network_igraph)$weight = as.numeric(pathway_network[,3])
	degs = degree(pathway_network_igraph, v = V(pathway_network_igraph))
	terminals = rep(0,length(V(pathway_network_igraph)))
	names(terminals) = V(pathway_network_igraph)$name
	#DEGs have positive prizes
	terminals[names(pathway_prizes)] = pathway_prizes
	#Steiner nodes have negative penalties as a function of their degree and alpha
	terminals[which(!names(terminals) %in% names(pathway_prizes))] = -degs[which(!names(terminals) %in% names(pathway_prizes))]^alpha
	#assign a high prize to the driver so it will be contained in the solution
	driver_prize = sum(pathway_prizes)*10
	terminals[driver_gene] = driver_prize
	res = PCSF(pathway_network_igraph, terminals = terminals, w = max(round(sum(pathway_prizes)/10),20), b = 1, mu = 0)
 	if(length(V(res)) == 0){ score = c(0,0)
	} else if(is_connected(res) & (driver_gene %in% V(res)$name)) {
		score = c(sum(V(res)$prize) - sum(E(res)$weight) -driver_prize -degs[driver_gene]^alpha,sum(V(res)$prize)-driver_prize)
	} else { score = c(-1,-1) }
	return(score)
}
