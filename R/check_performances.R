#' check_performances
#'
#' This function calculates mean precision, recall and F1 of the results provided by PRODIGY against a gold standard driver list
#' @param ranked_genes_lists A named list of ranked drivers for each patient.
#' @param snv_matrix A binary matrix with genes in rows and patients on columns. 1=mutation. All genes must be contained in the global PPI network.
#' @param gold_standard_drivers A list of known driver genes to be used as gold standard for validation.
#' @examples
#' data(gold_standard_drivers)
#' check_performances(results,snv_matrix,gold_standard_drivers)
#' @export
check_performances<-function(ranked_genes_lists,snv_matrix,gold_standard_drivers)
{
	libraries = c("ggplot2","cowplot")
	for(j in 1:length(libraries)){
	try({library(libraries[j],character.only=T)})
	}
	algorithm_colors = c("red")
	precision_matrices = list()
	recall_matrices = list()
	f1_matrices = list()	
	for(i in 1:length(ranked_genes_lists))
	{
		patient_snp = rownames(snv_matrix)[which(snv_matrix[,paste(names(ranked_genes_lists)[i],sep="")]==1)]
		if(length(intersect(gold_standard_drivers,patient_snp)) < 1) 
		{
			 print(paste("no known drivers with SNV mutations for patient",names(ranked_genes_lists)[i]))
			 next
		}
		if(length(intersect(gold_standard_drivers,ranked_genes_lists[[i]])) > 0)
		{
			precision_matrices[[names(ranked_genes_lists[i])]] = calculate_precision(ranked_genes_lists[[i]][1:min(20,length(ranked_genes_lists[[i]]))],gold_standard_drivers)[1:20]
			recall_matrices[[names(ranked_genes_lists[i])]] = calculate_recall(ranked_genes_lists[[i]][1:min(20,length(ranked_genes_lists[[i]]))],intersect(gold_standard_drivers,patient_snp))[1:20]
			curr_f1 = (2*precision_matrices[[names(ranked_genes_lists[i])]]*recall_matrices[[names(ranked_genes_lists[i])]])/(precision_matrices[[names(ranked_genes_lists[i])]]+recall_matrices[[names(ranked_genes_lists[i])]])
			curr_f1[is.nan(curr_f1)] = 0
			f1_matrices[[names(ranked_genes_lists[i])]] = curr_f1		
		}
	}
	precision_matrix = do.call(rbind,precision_matrices) 
	recall_matrix = do.call(rbind,recall_matrices)
	f1_matrix = do.call(rbind,f1_matrices)		
	return(list("precision"=precision_matrix,"recall"=recall_matrix,"f1"=f1_matrix))
	df_means = data.frame(x=seq(1,20,1))
	if(class(precision_matrix) == "matrix")
	{
		curr_precisions = apply(precision_matrix,2,function(x) mean(x))
	} else { curr_precisions = precision_matrix }
	df_means = cbind(df_means,curr_precisions,rep("PRODIGY",20))
	names(df_means)[c(2,3)] = c("value","Algorithm")
	precision_plot = ggplot(data=df_means,aes(x=x,y=value,colour=Algorithm))+  scale_colour_manual(values=algorithm_colors)+ geom_point(size=2) + geom_line(size=1) + labs(x="",y="Average precision") + ylim(c(0,ceiling(max(df_means[,"value"])*10)/10)) + theme(axis.title.y = element_text(size=26),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))
	
	df_means = data.frame(x=seq(1,20,1))
	if(class(recall_matrix) == "matrix")
	{
		curr_recalls = apply(recall_matrix,2,function(x) mean(x))
	} else { curr_recalls = recall_matrix }
	df_means = cbind(df_means,curr_recalls,rep("PRODIGY",20))
	names(df_means)[c(2,3)] = c("value","Algorithm")
	recall_plot = ggplot(data=df_means,aes(x=x,y=value,colour=Algorithm))+  scale_colour_manual(values=algorithm_colors)+ geom_point(size=2) + geom_line(size=1) + labs(x="",y="Average recall") + ylim(c(0,ceiling(max(df_means[,"value"])*10)/10)) + theme(axis.title.y = element_text(size=26),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))


	df_means = data.frame(x=seq(1,20,1))
	if(class(f1_matrix) == "matrix")
	{
		curr_f1 = apply(f1_matrix,2,function(x) mean(x))
	} else { curr_f1 = f1_matrix }
	df_means = cbind(df_means,curr_f1,rep("PRODIGY",20))
	names(df_means)[c(2,3)] = c("value","Algorithm")
	if(plot_figure)
	{
		f1_plot = ggplot(data=df_means,aes(x=x,y=value,colour=Algorithm))+  scale_colour_manual(values=algorithm_colors)+ geom_point(size=2) + geom_line(size=1) + labs(x="",y="Average F1") + ylim(c(0,ceiling(max(df_means[,"value"])*10)/10)) + theme(axis.title.y = element_text(size=26),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))
		plot_grid(title,precision_plot, recall_plot, f1_plot, nrow=4,ncol=1,rel_heights=c(0.1,1,1,1),rel_widths=c(1,1,1,1))
	}
}