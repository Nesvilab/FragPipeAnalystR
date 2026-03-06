#' Consensus clustering analysis
#'
#' Performs unsupervised consensus clustering to identify sample subtypes
#' using the ConsensusClusterPlus algorithm.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param additional_matrix Optional additional data matrix (e.g., RNA expression)
#'   to integrate with proteomics data. Can be a single matrix/data.frame or
#'   a list of matrices. Default is \code{NULL}.
#' @param max_cluster_num Numeric value for the maximum number of clusters to
#'   evaluate. Default is 10.
#' @param name Character string specifying the name for storing results in
#'   metadata. Default is "ConsensusClustering".
#' @param workdir Character string specifying the output directory for
#'   consensus clustering plots. Default is "ConsensusClusteringResult".
#' @param clusterAlg Character string specifying the clustering algorithm.
#'   Options include "pam" (default), "hc", "km", or others supported by
#'   ConsensusClusterPlus.
#' @param distance Character string specifying the distance metric.
#'   Options include "pearson" (default), "spearman", "euclidean", etc.
#' @param reps Numeric value for the number of resampling iterations.
#'   Default is 500.
#' @param pItem Numeric value between 0 and 1 for the proportion of items
#'   (samples) to sample in each iteration. Default is 0.8.
#' @param pFeature Numeric value between 0 and 1 for the proportion of
#'   features to sample in each iteration. Default is 1.
#' @param plot Character string specifying output format for plots.
#'   Options are "png" (default), "pdf", "ps", or NULL for no plots.
#' @param innerLinkage Character string for linkage method in inner clustering.
#'   Default is "average".
#' @param finalLinkage Character string for linkage method in final clustering.
#'   Default is "average".
#' @param writeTable Logical indicating whether to write consensus matrices
#'   to files. Default is \code{FALSE}.
#' @param weightsItem Optional numeric vector of weights for items.
#' @param weightsFeature Optional numeric vector of weights for features.
#' @param verbose Logical for verbose output. Default is \code{FALSE}.
#' @param corUse Character string for handling missing values in correlation.
#'   Default is "everything".
#' @param seed Numeric value for random seed. Default is \code{NULL} (uses 0).
#'
#' @return A \code{SummarizedExperiment} object with consensus clustering
#'   results stored in \code{metadata(se)[[name]]}.
#'
#' @details
#' Consensus clustering is a resampling-based method that provides robust
#' cluster assignments and helps determine the optimal number of clusters.
#' The function generates consensus matrices and diagnostic plots to evaluate
#' cluster stability.
#'
#' When \code{additional_matrix} is provided, the data are combined (rbind)
#' for integrated multi-omics clustering.
#'
#' @examples
#' \dontrun{
#' # Basic consensus clustering
#' se <- consensus_clustering_analysis(se, max_cluster_num = 6)
#'
#' # With custom parameters
#' se <- consensus_clustering_analysis(se, max_cluster_num = 8,
#'                                     clusterAlg = "hc",
#'                                     distance = "spearman",
#'                                     reps = 1000)
#'
#' # Multi-omics clustering
#' se <- consensus_clustering_analysis(se, additional_matrix = rna_matrix)
#' }
#'
#' @seealso \code{\link{snf_analysis}}
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom SummarizedExperiment assay metadata metadata<-
#'
#' @export
consensus_clustering_analysis <- function(se,
                                          additional_matrix=NULL, max_cluster_num=10, name="ConsensusClustering",
                                          workdir="ConsensusClusteringResult",
                                          clusterAlg="pam",
                                          distance="pearson",
                                          reps=500,
                                          pItem=0.8,
                                          pFeature=1,
                                          plot="png",
                                          innerLinkage="average",
                                          finalLinkage="average",
                                          writeTable=FALSE,
                                          weightsItem=NULL,
                                          weightsFeature=NULL,
                                          verbose=FALSE,
                                          corUse="everything",
                                          seed=NULL){
  if (is.null(seed)) {
    seed <- 0
  }

  if (is.null(additional_matrix)) {
    metadata(se)[[name]] <- list(result=ConsensusClusterPlus(
      assay(se),
      maxK=max_cluster_num,clusterAlg=clusterAlg,
      distance=distance,title=workdir,
      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
      innerLinkage=innerLinkage, finalLinkage=finalLinkage,
      writeTable=writeTable,weightsItem=weightsItem,weightsFeature=weightsFeature,
      verbose=verbose,corUse=corUse,seed=seed))
  } else {
    if(is.data.frame(additional_matrix) | is.matrix(additional_matrix)) {
      temp <- rbind(assay(se), additional_matrix)
    } else {
      temp <- assay(se)
      for(i in 1:length(additional_matrix)) {
        temp <- rbind(temp, additional_matrix[[i]])
      }
    }
    temp <- as.matrix(temp)
    metadata(se)[[name]] <- list(result=ConsensusClusterPlus(
      temp,
      maxK=max_cluster_num, clusterAlg=clusterAlg,
      distance=distance, title=workdir,
      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
      innerLinkage=innerLinkage, finalLinkage=finalLinkage,
      writeTable=writeTable, weightsItem=weightsItem,weightsFeature=weightsFeature,
      verbose=verbose,corUse=corUse,seed=seed))
  }
  return(se)
}


#' Similarity Network Fusion (SNF) analysis
#'
#' Performs multi-omics integration and clustering using Similarity Network
#' Fusion, which constructs sample similarity networks from each data type
#' and fuses them into a single unified network.
#'
#' @param se A \code{SummarizedExperiment} object containing proteomics data.
#' @param additional_matrix Additional data matrix (e.g., RNA expression) for
#'   multi-omics integration. Can be a single matrix/data.frame or a list of
#'   matrices. Default is \code{NULL}.
#' @param max_cluster_num Numeric value for the maximum number of clusters to
#'   evaluate. Default is 10.
#' @param name Character string specifying the name for storing results in
#'   metadata. Default is "SNFClustering".
#' @param workdir Character string specifying the output directory for
#'   SNF results and plots. Default is "SNF".
#' @param K Numeric value for the number of nearest neighbors to consider
#'   when constructing similarity graphs. Usually 10-30. Default is 20.
#' @param alpha Numeric value (0-1) controlling the kernel scaling parameter.
#'   Usually 0.3-0.8. Default is 0.5.
#' @param T Numeric value for the number of iterations in the diffusion process.
#'   Usually 10-20. Default is 10.
#' @param dist_metric Character string specifying the distance metric.
#'   Options are "euclidean" (default) or "pearson".
#' @param writeTable Logical indicating whether to write fusion graph and
#'   cluster results to files. Default is \code{FALSE}.
#' @param display_cluster Logical indicating whether to generate heatmap
#'   visualizations of the fused similarity matrix for each cluster number.
#'   Default is \code{TRUE}.
#'
#' @return A \code{SummarizedExperiment} object with SNF results stored in
#'   \code{metadata(se)[[name]]}, containing:
#'   \itemize{
#'     \item \code{dat_norm}: Normalized input data
#'     \item \code{dat_dist}: Distance matrices
#'     \item \code{dat_graph}: Affinity matrices for each data type
#'     \item \code{W}: Fused similarity network
#'     \item \code{result}: List of spectral clustering results for each k
#'   }
#'
#' @details
#' SNF is particularly useful for integrating multiple data types (e.g.,
#' proteomics and transcriptomics) where each data type contributes complementary
#' information about sample relationships. The method:
#' \enumerate{
#'   \item Normalizes each data type separately
#'   \item Constructs similarity networks for each data type
#'   \item Iteratively fuses networks using a diffusion process
#'   \item Performs spectral clustering on the fused network
#' }
#'
#' @examples
#' \dontrun{
#' # Single-omics SNF clustering
#' se <- snf_analysis(se, max_cluster_num = 5)
#'
#' # Multi-omics integration
#' se <- snf_analysis(se, additional_matrix = rna_matrix,
#'                    K = 15, alpha = 0.6)
#'
#' # Access clustering results for k=3
#' clusters_k3 <- metadata(se)$SNFClustering$result[[3]]$cluster
#' }
#'
#' @seealso \code{\link{consensus_clustering_analysis}}
#'
#' @importFrom SNFtool standardNormalization dist2 affinityMatrix SNF spectralClustering
#' @importFrom SummarizedExperiment assay colData metadata metadata<-
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grDevices png dev.off
#'
#' @export
snf_analysis = function(se,
                        additional_matrix=NULL,
                        max_cluster_num=10,
                        name="SNFClustering",
                        workdir="SNF",
                        K=20,
                        alpha=0.5,
                        T=10,
                        dist_metric="euclidean",
                        writeTable=F,
                        display_cluster=T) {
  min_cluster_num <- 2
  prefix <- "snf_graph"

  if (!file.exists(workdir)) {
    dir.create(workdir)
  }
  # K : number of neighbors, usually (10~30)
  # alpha : hyperparameter, usually (0.3~0.8)
  # T : Number of Iterations, usually (10~20)
  # C : number of clusters
  ## Input data is of size n x d_1, where n is the number of patients, d_1 is the number of genes.
  # Therefore, we need to transpose the data first
  if(is.data.frame(additional_matrix) | is.matrix(additional_matrix)) {
    dat <- list(t(assay(se)), t(additional_matrix))
  } else {
    dat <- list(t(assay(se)))
    for(i in 1: length(additional_matrix)) {
      dat <- append(dat, t(additional_matrix[[i]]))
    }
  }

  print("Perform SNF clustering analysis")
  ## If the data are all continuous values, we recommend the users to perform standard normalization before using SNF,
  ## though it is optional depending on the data the users want to use.
  message("*** do standard normalization")
  # use SNF with multiple views
  dat_norm = lapply(dat, function(x) SNFtool::standardNormalization(x))

  ## Calculate distance matrices (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
  ## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
  # if the data is discrete, we recommend the users to use ""
  # print(dat_norm)
  message("*** calculate distance matrices")
  if(dist_metric=="euclidean"){
    dat_dist = lapply(dat_norm, function(x) SNFtool::dist2(as.matrix(x), as.matrix(x))^(1/2))
  } else{
    dat_dist = lapply(dat_norm, function(x) cor(as.matrix(x), as.matrix(x), method = "pearson")) # correlation
  }

  ## next, construct similarity graphs
  message("*** construct similarity graphs")
  dat_graph = lapply(dat_dist, function(x) affinityMatrix(x, as.numeric(K), as.numeric(alpha)))
  # print(dat_graph)
  ## next, fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  message("*** fuse all the graphs")
  fusion_graph = SNF(dat_graph, K = K, t = T) # t: Number of iterations for the diffusion process.

  if (writeTable) {
    write.table(fusion_graph, file=file.path(workdir,paste0(prefix,"_fusion_graph.tsv")), row.names = TRUE, quote=F, sep = "\t")
  }

  dat_graph[["fusion_graph"]] = fusion_graph

  ## With this unified graph W (fusion_graph) of size n x n, you can do either spectral clustering or Kernel NMF.
  message("*** Working on spectral clustering")
  smp_lst_fused_lst = list()
  for(C in seq(min_cluster_num,max_cluster_num)){
    message(paste0("      clustering samples into ",C," groups"))
    # the final subtypes is determined by further clustering on fusion graph
    group <- spectralClustering(fusion_graph, C)

    if (writeTable) {
      cluster_result <- data.frame(cluster=group,
                                   Case_ID=colnames(fusion_graph))
      write.table(cluster_result, file.path(workdir, paste0(prefix, "_", C, "_clusters_result.tsv")),
                  row.names = F,sep="\t",quote=F)
    }

    normalize <- function(X) X/rowSums(X)
    ind <- sort(as.vector(group), index.return = TRUE)
    ind <- ind$ix
    W <- fusion_graph
    diag(W) <- median(as.vector(W))
    W <- normalize(W)
    W <- W + t(W)
    similarity <- W[ind, ind]

    names(group) <- colnames(se)
    smp_lst_fused_lst[[C]] <- list(matrix=similarity, cluster=group)

    if(display_cluster){
      png(file=file.path(workdir,paste0(prefix,"_SNF_", C, "_clusters.png")), height = 10, width = 10, units = "in", res = 300)
      draw(Heatmap(similarity, cluster_rows = F, cluster_columns = F), column_title=paste0("SNF k=", C))
      dev.off()
    }
  }

  metadata(se)[[name]] <- list(dat_norm=dat_norm,
                              dat_dist=dat_dist,
                              dat_graph=dat_graph,
                              W=fusion_graph,
                              result=smp_lst_fused_lst)
  return(se)
}
