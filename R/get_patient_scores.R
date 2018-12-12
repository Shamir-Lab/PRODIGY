get_patient_scores<-function(results_folder,diff_genes,patient_name,pathwayDB,pathwayDB_nodes)
{
	list_of_pathways = names(pathwayDB_nodes)
	pathways = list.files(results_folder)
	if(length(pathways) == 0)
	{
		return(NA)
	}
	scores = as.matrix(read.csv(paste(results_folder,pathways[1],sep=""),sep=" ",header=T,row.names=NULL))
	if(nrow(scores) < 2)
	{
		return(NA)
	}
	if(length(which(duplicated(scores[,1]))) > 0)
	{
		scores = scores[-which(duplicated(scores[,1])),]
	}
	rownames(scores) = scores[,1]
	scores = scores[,c(2,3)]
	scores[,1] = as.numeric(scores[,1])
	scores[,2] = as.numeric(scores[,2])
	total_scores = rep(0,nrow(scores))
	total_diff_scores = rep(0,nrow(scores))
	p_scores = matrix(ncol=length(rownames(scores)))
	colnames(p_scores) = rownames(scores)
	names(total_scores) = rownames(scores)
	names(total_diff_scores) = rownames(scores)	
	for(i in 1:length(pathways))
	{
		possible_error = tryCatch(
		{
				scores = as.matrix(read.csv(paste(results_folder,pathways[i],sep=""),sep=" ",header=T,row.names=NULL))
		},
		error=function(cond) {
				cond
		})
		if(inherits(possible_error, "error")){ next }
		if(length(which(duplicated(scores[,1]))) > 0)
		{
			scores = scores[-which(duplicated(scores[,1])),]
		}
		scores_2 = matrix(ncol=2,nrow = nrow(scores))			
		rownames(scores_2) = scores[,1]
		scores_2[,1] = as.numeric(scores[,2])
		scores_2[,2] = as.numeric(scores[,3])
		scores = scores_2
		pathway_name = substr(pathways[i],1,regexpr(pattern = "_",pathways[i])[1]-1)
		if(pathwayDB == "nci") 
		{ 
			if(pathway_name == "TCR signaling in naive CD4+ T cells"){ pathway_name =  "TCR signaling in na?ve CD4+ T cells"}
			if(pathway_name == "TCR signaling in naive CD8+ T cells"){ pathway_name =  "TCR signaling in na?ve CD8+ T cells"}
			if(pathway_name == "Downstream signaling in na?ve CD8+ T cells"){ pathway_name =  "Downstream signaling in na?ve CD8+ T cells"}				
		}
		pathway_name = gsub("%","/",pathway_name)
		pathway_name = gsub(";",":",pathway_name) 					
		if(!pathway_name %in% list_of_pathways){ next }
		seed = intersect(names(diff_genes),pathwayDB_nodes[[pathway_name]])
		max_score = sum(diff_genes[seed])
		if(is.na(max_score))
		{
			next
		}
		scores[,1] = round(as.numeric(scores[,1]),digits=3)
		scores[,1] = scores[,1]/max_score
		scores[,2] = scores[,2]/max_score
		p_scores = rbind(p_scores,scores[,1])
		rownames(p_scores)[nrow(p_scores)] = substr(pathways[i],1,regexpr(pattern = "_scores",pathways[i])[1]-1)				
	}
	return(p_scores[-1,])
}
