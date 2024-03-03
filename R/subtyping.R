#' @param se input SummarizedExperiment object
#' @param additional_matrix additional matrix like RNA expression
#' @param max_cluster_num
#' @param name
#' @param workdir
#' @param clusterAlg
#' @param distance
#' @param reps
#' @param pItem
#' @param pFeature
#' @param plot description
#' @param innerLinkage
#' @param finalLinkage
#' @param writeTable
#' @param weightsItem
#' @param weightsFeature
#' @param verbose
#' @param corUse
#' @param seed
#' @return SummarizedExperiment object with consensus clustering analysis result added into metadata
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


#' @param se input SummarizedExperiment object
#' @param additional_matrix additional matrix like RNA expression
#' @param max_cluster_num
#' @param K
#' @param alpha
#' @param dist_metric
#' @param cluster_method
#' @param display_cluster
#' @param workdir
#'
#' @return SummarizedExperiment object with SNF analysis result added into metadata
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
