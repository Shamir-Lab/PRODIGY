#' analyze_PRODIGY_results
#'
#' This function analyzes PRODIGY's results (influence score matrices) to rank individual driver genes for a group of patients.
#' @param all_patients_scores Either a list of influence scores matrices or a single matrix received by the PRODIGY algorithm.
#' @return A list of ranked genes per sample.
#' @examples
#' results = analyze_PRODIGY_results(all_patients_scores)
#' @export
analyze_PRODIGY_results<-function(all_patients_scores)
{
	libraries = c("DESeq2","igraph","plyr","biomaRt","MASS","mixtools")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	if(class(all_patients_scores) == "matrix"){ all_patients_scores = list(all_patients_scores) }
	for(i in 1:length(all_patients_scores))
	{	
		all_patients_scores[[i]][which(is.na(all_patients_scores[[i]]))] = 0
		all_patients_scores[[i]][all_patients_scores[[i]] < 0] = 0
	}
	Prodigy_rankings = list()
	for(i in 1:length(all_patients_scores))
	{
		if(is.null(all_patients_scores[[i]])){Prodigy_rankings[[i]] = c(); next}
		pathways_to_take = c()
		#single pathway
		if(is.null(nrow(all_patients_scores[[i]])))
		{
			ranking = sort(all_patients_scores[[i]],decreasing=T)
		#less than 4 mutations, no pathway filtering
		} else if(ncol(all_patients_scores[[i]]) < 4) {
			ranking = sort(apply(all_patients_scores[[i]],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)	
		#filter pathways
		} else {
			pathways_to_take = names(which(apply(all_patients_scores[[i]],1,function(x) length(which(x > 0))) < ncol(all_patients_scores[[i]])/2))
			if(length(pathways_to_take) < 2)
			{
				ranking = sort(all_patients_scores[[i]][pathways_to_take,],decreasing=T)
				#if all pathways are filtered, abort pathway filtering
				if(all(ranking==0))
				{
					pathways_to_take = rownames(all_patients_scores[[i]])
					ranking = sort(apply(all_patients_scores[[i]][pathways_to_take,],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				}
			} else {
				ranking = sort(apply(all_patients_scores[[i]][pathways_to_take,],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				if(all(ranking==0))
				{
					pathways_to_take = rownames(all_patients_scores[[i]])
					ranking = sort(apply(all_patients_scores[[i]][pathways_to_take,],2,function(x) sum(sort(x,decreasing=T))),decreasing=T)
				} else {
					pathways_to_take = pathways_to_take[which(apply(all_patients_scores[[i]][pathways_to_take,],1,sum)!=0)]
				}
			}
		}
		ranking = ranking[ranking > 0]
		if(length(ranking) == 0) {Prodigy_rankings[[i]] = c(); next}
		# single pathway was used, no bimodel distribution 
		if(length(pathways_to_take) < 2){Prodigy_rankings[[i]] = names(ranking);next}
		# check bimodel distribution
		possible_error = tryCatch(
		{
				bimodel_dist = normalmixEM(ranking,k=2)
				unimodel_dist = MASS::fitdistr(ranking,"normal")
		},
		error=function(cond) {
				cond
		})
		if(inherits(possible_error, "error")){ Prodigy_rankings[[i]] = names(ranking)
		} else {
			# check which distribution is more likely
			if(bimodel_dist$loglik > unimodel_dist$loglik)
			{
					Prodigy_rankings[[i]] = names(ranking[bimodel_dist$posterior[,which(bimodel_dist$mu == max(bimodel_dist$mu))] >0.5])
			} else {
					Prodigy_rankings[[i]] = names(ranking)
			}
		}
	}
	names(Prodigy_rankings) = names(all_patients_scores)
	return(Prodigy_rankings)
}