#' @export
get_diff_expressed_genes<-function(expression_matrix,sample,sample_origins,beta=2,gamma=0.05)
{
	if(length(which(sample_origins == "normal")) < 1)
	{	
			print("no normal samples, cannot perform differential expression analysis. aborting")
			return()
	}
	df = data.frame(expression_matrix[,sample])
	healthy_tissues = which(sample_origins=="normal")
	df = cbind(df,expression_matrix[,healthy_tissues])
	colnames(df)[1] = paste(sample)
	coldata = data.frame(rep("healthy",ncol(df)))
	colnames(coldata) = "condition"
	rownames(coldata) = colnames(df)
	coldata$condition = factor(coldata$condition, levels = c("healthy","tumor"))
	coldata[1,"condition"] = "tumor"
	#countData is a data frame with genes on rows and patients on columns
	#coldata is a data frame with patients on the rows and "tumor"/"healthy" on the first row
	diff_expression_count_matrix = DESeqDataSetFromMatrix(countData = round(df),
                              colData = coldata,
                              design = ~ condition)
	diff_expression_count_matrix = diff_expression_count_matrix[rowSums(counts(diff_expression_count_matrix)) > 1,]
	diff_expression_count_matrix$condition = factor(diff_expression_count_matrix$condition, levels = c("tumor","healthy"))
	diff_expression_count_matrix = DESeq(diff_expression_count_matrix,betaPrior=TRUE)
	res = results(diff_expression_count_matrix,alpha=gamma,lfcThreshold=1,altHypothesis="greaterAbs")
	res = res[!is.na(res$log2FoldChange),]
	res = res[!is.na(res$pvalue),]
	diff_genes = res[res$pvalue<gamma & abs(res$log2FoldChange) > beta,2]
	names(diff_genes) = rownames(res)[res$pvalue<gamma & abs(res$log2FoldChange) > beta]
	return(abs(diff_genes))
}