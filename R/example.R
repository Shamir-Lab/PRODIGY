#' Example run of PRODIGY
#' download the package and change working directory to PRODIGY's root folder
#' install("PRODIGY")
#' library(PRODIGY)
#' # Load SNP+expression data from TCGA
#' data(COAD_SNP)
#' data(COAD_Expression)
#' # Load STRING network data 
#' data(STRING_network)
#' network = STRING_network
#' # Take samples for which SNP and expression is available 
#' samples = intersect(colnames(expression_matrix),colnames(snp_matrix))[1:5]
#' # Get differentially expressed genes (DEGs) for all samples
#' expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
#' DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL)
#' # Identify sample origins (tumor or normal)
#' sample_origins = rep("tumor",ncol(expression_matrix))
#' sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
#' # Run PRODIGY
#' all_patients_scores = PRODIGY_cohort(snp_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,
#' 			pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins)
#' # Get driver gene rankings for all samples 
#' results = analyze_PRODIGY_results(all_patients_scores) 