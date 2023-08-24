#

#Space Time
#Author: DD
#E-mail:fngseng12345@163.com


#' Split data matrix into smaller sub-matrices ('chunks')
#'
#' @param   matrix      Input data matrix
#' @param   chunk.size  How many cells to include in each sub-matrix
#'
#' @return  A list of sub-matrices, each with size {n_features x chunk_size}
#'
split_data.matrix <- function(matrix, chunk.size=1000) {
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1

  split.data <- list()
  min <- 1
  for (i in seq_len(nchunks)) {
    if (i == nchunks-1) {  #make last two chunks of equal size
      left <- ncols-(i-1)*chunk.size
      max <- min+round(left/2)-1
    } else {
      max <- min(i*chunk.size, ncols)
    }
    split.data[[i]] <- matrix[,min:max]
    min <- max+1    #for next chunk
  }
  return(split.data)
}

#' run_VISION_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_VISION_pipline <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="counts",BPPARAM=NULL,ncores=1){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }

  matrix <- Seurat::GetAssayData(object=rds, slot=slot, assay=assay)

  # use counts,normalize matrix
  if (slot=="counts"){
    n.umi <- colSums(matrix)
    center.umi <- median(n.umi)
    matrix_list <-  split_data.matrix(matrix,chunk.size=chunk.size)

    scaled_counts_list <- lapply(matrix_list,function(x) {
      x <- as.matrix(x)
      n.umi <- colSums(x)
      scaled_counts <- t(t(x) / n.umi) * center.umi
      return(scaled_counts)
    }
    )
  }

  if (slot=="data"){
    scaled_counts_list <- split_data.matrix(matrix,chunk.size=chunk.size)
  }

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }

  meta.list <- BiocParallel::bplapply(
    X = scaled_counts_list,
    BPPARAM =  BPPARAM,
    FUN = function(x) {
      set.seed(123)
      vis <- VISION::Vision(x, signatures = signatures,min_signature_genes = 0.001)
      options(mc.cores = ncores)
      vis <- VISION::analyze(vis)
      signature_exp <- data.frame(vis@SigScores)
      return(list(signature_exp=signature_exp))
    }
  )
  meta.merge <- lapply(meta.list,function(x) cbind(x[["signature_exp"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' run_AUCell_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_AUCell_pipline <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="counts",BPPARAM=NULL,ncores=1){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }

  matrix <- Seurat::GetAssayData(object=rds, slot=slot, assay=assay)

  scaled_counts_list <- split_data.matrix(matrix,chunk.size=chunk.size)

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }

  geneSets <- GSEABase::getGmt(signatures)

  meta.list <- BiocParallel::bplapply(
    X = scaled_counts_list,
    BPPARAM =  BPPARAM,
    FUN = function(x) {
      set.seed(123)
      cells_rankings <- AUCell::AUCell_buildRankings(as.matrix(x), nCores=ncores, plotStats=F)
      cells_AUC <- AUCell::AUCell_calcAUC(geneSets, cells_rankings)
      signature_exp <- data.frame(t(AUCell::getAUC(cells_AUC)))
      return(list(signature_exp=signature_exp))
    }
  )
  meta.merge <- lapply(meta.list,function(x) cbind(x[["signature_exp"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' run_GSVA_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param method gsva method,'gsva', 'ssgsea', 'zscore', 'plage'
#' @param kcdf  "Gaussian", "Poisson", "none"   Character string denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples when ‘method="gsva"’.  By default, ‘kcdf="Gaussian"’ which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.  When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to ‘kcdf="Poisson"’.
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_GSVA_pipline <- function(rds,signatures,method="ssgsea",chunk.size=1000,assay=NULL,slot="counts",kcdf="Poisson",BPPARAM=NULL,ncores=1){
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }

  matrix <- Seurat::GetAssayData(object=rds, slot=slot, assay=assay)

  scaled_counts_list <- split_data.matrix(matrix,chunk.size=chunk.size)

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }

  geneSets <- GSEABase::getGmt(signatures)

  meta.list <- BiocParallel::bplapply(
    X = scaled_counts_list,
    BPPARAM =  BPPARAM,
    FUN = function(x) {
      set.seed(123)
      gsva_es <- GSVA::gsva(as.matrix(x), geneSets, method=c(method), tau=switch(method, gsva=1, ssgsea=0.25, NA),kcdf=kcdf, parallel.sz=ncores) #
      signature_exp<-data.frame(t(gsva_es))
      return(list(signature_exp=signature_exp))
    }
  )
  meta.merge <- lapply(meta.list,function(x) cbind(x[["signature_exp"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' run_UCell_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_UCell_pipline <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="data",BPPARAM=NULL,ncores=1){
  geneSets <- GSEABase::getGmt(signatures)
  markers <- GSEABase::geneIds(geneSets)
  rds1 <- UCell::AddModuleScore_UCell(rds, features = markers,assay=assay,slot=slot,BPPARAM=BPPARAM,ncores=ncores,chunk.size=chunk.size)
  metadata <- rds1@meta.data

  pathway <- c()
  for (i in names(markers)){
    pathway <- c(pathway,colnames(metadata)[grep(i,colnames(metadata))])
  }
  meta.merge <- Seurat::FetchData(rds1,vars=pathway)
  colnames(meta.merge) <- gsub("_UCell","",colnames(meta.merge))
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' run_AddModuleScore_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_AddModuleScore_pipline <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="data",BPPARAM=NULL,ncores=1){
  geneSets <- GSEABase::getGmt(signatures)
  markers <- GSEABase::geneIds(geneSets)
  # names(markers) <- gsub("-","_",names(markers))
  # rds1 <-  AddModuleScore(object = rds,features = markers,name =names(markers))
  pathname <- "pathwayname"
  rds1 <-  Seurat::AddModuleScore(object = rds,features = markers,name =pathname)
  metadata <- rds1@meta.data

  # pathway <- c()
  # for (i in names(markers)){
  # pathway <- c(pathway,colnames(metadata)[grep(i,colnames(metadata))])
  # }
  pathway <- colnames(rds1@meta.data)[grep(pathname,colnames(rds1@meta.data))]
  meta.merge <- Seurat::FetchData(rds1,vars=pathway)
  colnames(meta.merge) <-names(markers)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' Addsingscore
#'
#' @param matrix  single cell expression matrix
#' @param signatures gmt file path
#'
#' @return  matrix
#' @export
#'
#' @examples
#'
#'
Addsingscore <- function(matrix,signatures){
  geneSets <- GSEABase::getGmt(signatures)
  markers <- GSEABase::geneIds(geneSets)
  h.gsets.list <- markers %>% purrr::compact()

  singscore.rank <- singscore::rankGenes(as.data.frame(matrix))
  # calculate separately
  singscore.scores <- list()
  for (i in seq_along(h.gsets.list)){
    if (any(stringr::str_detect(h.gsets.list[[i]], pattern = "\\+$|-$"))) {
      h.gsets.list.positive <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)\\+")[,2] %>% purrr::discard(is.na)
      h.gsets.list.negative <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)-")[,2] %>% purrr::discard(is.na)
      if (length(h.gsets.list.positive)==0) {
        singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,upSet = h.gsets.list.negative, centerScore = F)
      }
      if (length(h.gsets.list.negative)==0) {
        singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,upSet = h.gsets.list.positive, centerScore = F)
      }
      if ((length(h.gsets.list.positive)!=0)&(length(h.gsets.list.negative)!=0)) {
        singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,upSet = h.gsets.list.positive,downSet = h.gsets.list.negative, centerScore = F)
      }
    }else{
      singscore.scores[[i]] <- singscore::simpleScore(singscore.rank, upSet = h.gsets.list[[i]], centerScore = F)
    }
    TotalScore <- NULL
    singscore.scores[[i]] <- singscore.scores[[i]] %>%
      dplyr::select(TotalScore) %>%
      magrittr::set_colnames(names(h.gsets.list)[i])}
  names(singscore.scores) <- names(h.gsets.list)
  singscore.scores <- do.call(cbind, singscore.scores)
  return(singscore.scores)
}

#' run_singscore_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

run_singscore_pipline <- function(rds,signatures,chunk.size=1000,assay=NULL,slot="data",BPPARAM=NULL,ncores=1){

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }

  matrix <- Seurat::GetAssayData(object=rds, slot=slot, assay=assay)

  scaled_counts_list <- split_data.matrix(matrix,chunk.size=chunk.size)

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }

  meta.list <- BiocParallel::bplapply(
    X = scaled_counts_list,
    BPPARAM =  BPPARAM,
    FUN = function(x) {
      set.seed(123)
      signature_exp <- Addsingscore(x, signatures)
      return(list(signature_exp=signature_exp))
    }
  )
  meta.merge <- lapply(meta.list,function(x) cbind(x[["signature_exp"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' RankCalculation: Calculating Mean Ranks for signature genes across each cell
#'
#' @param x  single cell expression matrix
#' @param genes
#'
#' @return  matrix
#' @export
#'
#' @examples
#'
#'

RankCalculation <- function(x,genes){

  subdata = x[x!=0]                                                                      ### Removing Dropouts from single cell
  DataRanksUpdated=rank(subdata)                                                         ### Calculating ranks of each signature gene per cell
  DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]        ### Shortling rank vector for signature genes
  CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 )     ### Calculating Mean of ranks for signature genes
  FinalRawRank = CumSum/length(subdata)                                                  ### Normalizing Means by total coverage
  return(FinalRawRank)
}

#' ORCalculation: Calculating enrichment of signature genes across each cell 	(using odds ratio)
#'
#' @param data  single cell expression matrix
#' @param genes
#'
#' @return  matrix
#' @export
#'
#' @examples
#'
#'

ORCalculation <- function(data,genes){
  GE = data[which(rownames(data) %in% genes),]                                          ### Subsetting data for signature genes
  NGE = data[-which(rownames(data) %in% genes),]                                        ### Subsetting data for non-signature genes
  SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))                                 ### Calculating Number of expressed Signature Genes per cell
  NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))                               ### Calculating Number of expressed Non-Signature Genes per cell
  SigGenesNE = nrow(GE) - SigGenesExp                                                   ### Calculating Number of Not expressed Signature Genes per cell
  SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)									  ### Replacing Zero's with 1
  NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)                                ### Replacing Zero's with 1
  NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)                               ### Calculating Number of Not expressed Non-Signature Genes per cell
  NSigGenesNE = NSigGenesNE - SigGenesNE
  OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)                         ### Calculating Enrichment (Odds Ratio)
  return(OR)
}

#' LikelihoodCalculation: Calculating enrichment of signature genes across each cell (using Likelihood ratio)
#'
#' @param data  single cell expression matrix
#' @param genes
#'
#' @return  matrix
#' @export
#'
#' @examples
#'
#'

LikelihoodCalculation <- function(data,genes){
  GE = data[which(rownames(data) %in% genes),]
  NGE = data[-which(rownames(data) %in% genes),]
  SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
  NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
  SigGenesNE = nrow(GE) - SigGenesExp
  SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)
  NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
  NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
  NSigGenesNE = NSigGenesNE - SigGenesNE
  LR1 = SigGenesExp*(NSigGenesExp + NSigGenesNE)
  LR2 = NSigGenesExp * (SigGenesExp + SigGenesNE)
  LR = LR1/LR2
  return(LR)
}

#' NormalizationJAS: Scalar [0,1] Normalization of Means and Enrichment across set of cells
#'
#' @param JAS_Scores
#'
#' @return  JAS_Scores
#' @export
#'
#' @examples
#'
#'
NormalizationJAS <- function(JAS_Scores)
{
  JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
  return(JAS_Scores)
}

#' JASMINE: Signature Scoring via JASMINE mergining Means and Enrichment
#'
#' @param data
#' @param genes
#' @param method oddsratio or likelihood
#'
#' @return  JAS_Scores
#' @export
#'
#' @examples
#'
#'
JASMINE <- function(data,genes,method)
{
  idx = match(genes,rownames(data))
  idx = idx[!is.na(idx)]
  if(length(idx)> 1){
    RM = apply(data,2,function(x) RankCalculation(x,genes))                              ### Mean RankCalculation for single cell data matrix
    RM = NormalizationJAS(RM)                                                            ### Normalizing Mean Ranks

    if(method == "oddsratio"){
      OR = ORCalculation(data,genes)			                                             ### Signature Enrichment Calculation for single cell data matrix (OR)
      OR = NormalizationJAS(OR)															 ### Normalizing Enrichment Scores (OR)
      JAS_Scores = (RM + OR)/2
    }else if(method == "likelihood"){

      LR = LikelihoodCalculation(data,genes)			                                     ### Signature Enrichment Calculation for single cell data matrix  (LR)
      LR = NormalizationJAS(LR)															 ### Normalizing Enrichment Scores (LR)
      JAS_Scores = (RM + LR)/2
    }
    FinalScores = data.frame(names(RM),JAS_Scores)                                       ### JASMINE scores
    colnames(FinalScores)[1]='SampleID'
    return(FinalScores)
  }
}

#' AddModuleScore_JAS
#'
#' @param matrix
#' @param signatures
#' @param method oddsratio or likelihood
#'
#' @return  matrix
#' @export
#'
#' @examples
#'
#'
AddModuleScore_JAS  <- function(matrix,signatures,method){
  geneSets <- clusterProfiler::read.gmt(signatures)
  pathway <- levels(geneSets$term)
  geneSets$term <- as.character(geneSets$term)
  geneSets$gene <- as.character(geneSets$gene)
  i <- pathway[1]
  gene <- geneSets[which(geneSets$term == i),]$gene
  JAS_result_tmp  <-   JASMINE(matrix,gene,method = method)
  JAS_result_tmp <- JAS_result_tmp[,2,drop=FALSE]
  colnames(JAS_result_tmp) <- i
  JAS_result <- JAS_result_tmp
  if (length(pathway)>1){
    for (i in pathway[2:length(pathway)]){
      gene <- geneSets[which(geneSets$term == i),]$gene
      print(i)
      print(gene)
      JAS_result_tmp  <-   JASMINE(matrix,gene,method = method)
      JAS_result_tmp <- JAS_result_tmp[,2,drop=FALSE]
      colnames(JAS_result_tmp) <- i
      JAS_result <- cbind(JAS_result,JAS_result_tmp)
    }
  }
  return(JAS_result)
}

#' run_JAS_pipline
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object,(if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation. If \code{BPPARAM = NULL}, the function uses,\code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param method oddsratio or likelihood
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'
run_JAS_pipline <- function(rds,slot="counts",assay=NULL,chunk.size=1000,signatures,BPPARAM=NULL,ncores=1,method="likelihood"){


  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(rds)
  }

  matrix <- Seurat::GetAssayData(object=rds, slot=slot, assay=assay)

  scaled_counts_list <- split_data.matrix(matrix,chunk.size=chunk.size)

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }

  meta.list <- BiocParallel::bplapply(
    X = scaled_counts_list,
    BPPARAM =  BPPARAM,
    FUN = function(x) {
      set.seed(123)
      signature_exp <- AddModuleScore_JAS(x,signatures,method)
      return(list(signature_exp=signature_exp))
    }
  )
  meta.merge <- lapply(meta.list,function(x) cbind(x[["signature_exp"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  rds <- Seurat::AddMetaData(rds, as.data.frame(meta.merge))
  return(rds)
}

#' scgmt
#'
#' @param rds      seurat object
#' @param chunk.size  How many cells to include in each sub-matrix
#' @param signatures gmt file path
#' @param assay Pull out data from this assay of the Seurat object
#' @param slot Pull out data from this slot of the Seurat object,default is counts
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param ncores ncores Number of processors to parallelize computation.
#' @param kcdf
#' @param method "AUCell","UCell","AddModuleScore","gsva","ssgsea","zscore","plage","VISION","JAS_likelihood","JAS_oddsratio","singscore"
#'
#' @return  seurat object
#' @export
#'
#' @examples
#'
#'

scgmt <- function(rds,slot=NULL,assay=NULL,chunk.size=1000,signatures,BPPARAM=NULL,kcdf =NULL,ncores=1,method){
  if (method %in% c("AUCell","UCell","AddModuleScore","gsva","ssgsea","zscore","plage","VISION","JAS_likelihood","JAS_oddsratio","singscore")){

    print(paste0("------------------run ",method,"------------------"))
    timestart <- Sys.time()
    if (method=="AUCell"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_AUCell_pipline(rds=rds,chunk.size=chunk.size,assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="UCell"){
      if (is.null(slot)) {
        slot <- "data"
      }
      rds1 <- run_UCell_pipline(rds=rds,chunk.size=chunk.size,assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="AddModuleScore"){
      if (is.null(slot)) {
        slot <- "data"
      }
      rds1 <- run_AddModuleScore_pipline(rds=rds,chunk.size=chunk.size,assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method %in% c('gsva', 'ssgsea', 'zscore', 'plage')){
      if (is.null(slot)) {
        slot <- "counts"
      }
      if (is.null(kcdf)) {
        if (method %in% c('gsva', 'ssgsea')){
          kcdf <- "Poisson"
        }
        if (method %in% c('zscore', 'plage')){
          kcdf <- "Gaussian"
        }
      } else {
        kcdf <- kcdf
      }

      rds1 <- run_GSVA_pipline(rds=rds,chunk.size=chunk.size,assay=assay,method=method,slot=slot,signatures=signatures,kcdf=kcdf,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="VISION"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_VISION_pipline(rds=rds,chunk.size=chunk.size,assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="JAS_likelihood"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_JAS_pipline(rds=rds,chunk.size=chunk.size,method="likelihood",assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="JAS_oddsratio"){
      if (is.null(slot)) {
        slot <- "counts"
      }
      rds1 <- run_JAS_pipline(rds=rds,chunk.size=chunk.size,method="oddsratio",assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    if (method=="singscore"){
      if (is.null(slot)) {
        slot <- "data"
      }
      rds1 <- run_singscore_pipline(rds=rds,chunk.size=chunk.size,assay=assay,slot=slot,signatures=signatures,BPPARAM=BPPARAM,ncores=ncores)
    }
    timeend <- Sys.time()
    time_use <- round(difftime(timeend,timestart,units = "secs"),2)
    print(paste0("time use :",time_use," s"))
    print(paste0("------------------run ",method," done!------------------"))

    return(rds1)
  } else {
    stop(paste0("please check method !"))
  }
}

#' scgmt_line_plot
#'
#' @param object      seurat object
#' @param signatures gmt file path
#' @param group.by
#' @param x.lab.title
#' @param y.lab.title
#' @param face
#' @param line.size
#'
#' @return  ggplot2 object
#' @export
#'
#' @examples
#'
#'

scgmt_line_plot <- function(object,signatures,group.by=NULL,x.lab.title=NULL,y.lab.title=NULL,face=NULL,line.size=1.2){

  if (methods::is(object)=="Seurat"){
    meta <- object@meta.data
  } else if (class(object)=="data.frame"){
    meta <- object
  }

  if (is.null(x.lab.title)){
    x.lab.title <- paste0("The score of ",signatures)
  } else {
    x.lab.title <- x.lab.title
  }

  if (is.null(y.lab.title)){
    y.lab.title <- "Cumulative Distribution Probability"
  } else {
    y.lab.title <- y.lab.title
  }

  p <- ggplot2::ggplot(data=meta,aes(x=meta[,signatures],color=meta[,group.by])) +
    stat_ecdf(size=line.size) +
    theme_classic() + labs(x=x.lab.title,y= y.lab.title,color=group.by) +
    theme(plot.title=element_text(hjust=0.5),
          axis.title.x =element_text (face = face,size=10),
          axis.title.y =element_text (face = face,size=10),
          legend.title = element_text(size = 10),
          legend.text = element_text(face = face,size = 10),
          axis.text.y = element_text (face = face,size=10),
          axis.text.x = element_text (face = face,size=10),
          axis.text = element_text (color = "black"))
  return(p)
}

#' scgmt_merge_line_plot
#'
#' @param object      seurat object
#' @param signatures gmt file path
#' @param x.lab.title
#' @param y.lab.title
#' @param face
#' @param line.size
#'
#' @return  ggplot2 object
#' @export
#'
#' @examples
#'
#'

scgmt_merge_line_plot <- function(object,signatures,x.lab.title=NULL,y.lab.title=NULL,face=NULL,line.size=1.2){

  meta <- object@meta.data
  signatures_merge <- intersect(signatures,colnames(meta))
  ll <- setdiff(signatures,signatures_merge)
  if (length(ll) > 0.5){
    mm <- paste0("Some pathway is not in data: ",ll)
    warning(mm, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  meta_sub <- meta[,signatures_merge,drop=F]
  data <- reshape2::melt(meta_sub)
  colnames(data) <- c("pathway","score")
  p <- scgmt_line_plot(object=data,signatures="score",group.by="pathway",x.lab.title="The score of pathway",face=face,line.size=line.size)
  return(p)
}

#' scgmt_heatmap_plot
#'
#' @param object      seurat object
#' @param signatures gmt file path
#' @param group.by
#' @param cluster_cols
#' @param cluster_rows
#' @param scale
#' @param group.order
#' @param pathway.order
#' @param heat_col
#' @param border_color
#'
#' @return  pheatmap object
#' @export
#'
#' @examples
#'
#'

scgmt_heatmap_plot <- function(object,signatures,group.by="ident",cluster_cols=FALSE,cluster_rows=FALSE,scale="row",group.order=NULL,pathway.order=NULL,heat_col=NULL,border_color="white",...){
  object$cluster_cell <- object@active.ident
  meta <- object@meta.data
  signatures_merge <- intersect(signatures,colnames(meta))
  ll <- setdiff(signatures,signatures_merge)
  if (length(ll) > 0.5){
    mm <- paste0("Some pathway is not in data: ",ll)
    warning(mm, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  if (group.by=="ident"){
    group.by <- "cluster_cell"
  } else {
    group.by <- group.by
  }
  subclomn <- c(group.by,signatures_merge)
  subclomn_new <- c("group_ident",signatures_merge)
  meta_sub <- meta[,subclomn,drop=F]
  colnames(meta_sub) <- subclomn_new
  meta_ave <- aggregate(.~group_ident,meta_sub,mean)
  rownames(meta_ave) <-  meta_ave$group_ident
  meta_ave <- meta_ave[,signatures_merge,drop=F]
  meta_ave <- t(meta_ave)

  if (!is.null(group.order)) {
    meta_ave <- meta_ave[,group.order,drop=F]
  }

  if (!is.null(pathway.order)) {
    meta_ave <- meta_ave[pathway.order,,drop=F]
  }

  if(is.null(heat_col)){
    heatmap_col <- c("#0099CC", "white", "#CC0033")
  }else{
    heatmap_col <- heat_col
  }
  col_fun <- colorRampPalette(heatmap_col)(100)
  heatmapfont <- 20
  p <- pheatmap::pheatmap(meta_ave,color = col_fun,cluster_cols=cluster_cols,cluster_rows=cluster_rows,scale = scale,cellwidth =heatmapfont+2,cellheight=heatmapfont+2,border_color=border_color,fontsize=heatmapfont,...)
  return(p)
}

#' scgmt_scatter_plot
#'
#' @param rds      seurat object
#' @param signatures.x
#' @param signatures.y
#' @param group.by
#' @param cols
#' @param pt.size
#' @param guideline
#'
#' @return  ggplot2 object
#' @export
#'
#' @examples
#'
#'

scgmt_scatter_plot <- function(rds,signature.x,signature.y,cols=NULL,group.by="ident",pt.size=1,guideline=TRUE){

  rds$ident <- rds@active.ident
  meta <- rds@meta.data
  if (!signature.x %in% colnames(meta)){stop(paste0("signature.x: ",signature.x," is not in data,please check !"))}
  if (!signature.y %in% colnames(meta)){stop(paste0("signature.y: ",signature.y," is not in data,please check !"))}
  if (!group.by %in% colnames(meta)){stop(paste0("group.by: ",group.by," is not in data,please check !"))}
  data <- meta[,c(signature.x,signature.y,group.by)]
  colnames(data) <- c("signature.x","signature.y","group")
  center.min <- round(min(min(data$signature.x),min(data$signature.y)), digits = 1) -0.1
  if (is.null(cols)){
    p  <-   ggplot2::ggplot(data,aes(x=signature.x,y=signature.y,color=group)) +
      geom_point(size=pt.size)
  } else {
    p  <-   ggplot2::ggplot(data,aes(x=signature.x,y=signature.y,color=group)) +
      geom_point(size=pt.size) +
      scale_color_manual(values=cols)
  }
  p  <-   p +
    theme_classic() +
    expand_limits(x = center.min, y = center.min) +
    # scale_x_continuous(limits = c(-0.2,1),
    # breaks = seq(-0.2,1,by=0.2))+
    # scale_y_continuous(limits = c(-0.2,1),
    # breaks = seq(-0.2,1,by=0.2)) +
    labs(x=signature.x,y=signature.y,color="") +
    theme(axis.text=element_text(color="black"),
          axis.text.x=element_text(size=8,face="bold"),
          axis.text.y=element_text(size=8,face="bold"),
          plot.title = element_text(size=14,hjust = 0.5,color="black",face="bold"),
          plot.subtitle=element_text(hjust = 0.5,color="black",face="bold"),
          axis.title=element_text(size=12,face="bold") )
  if (guideline) {
    p <- p + geom_abline(slope = 1,intercept = 0,lty="dashed")
  }

  return(p)
}

#' scgmt_ridges_plot
#'
#' @param rds      seurat object
#' @param signature
#' @param group.by
#'
#' @return  ggplot2 object
#' @export
#'
#' @examples
#'
#'

scgmt_ridges_plot <- function(rds,signature,group.by="ident"){

  rds$ident <- rds@active.ident
  meta <- rds@meta.data
  if (!signature %in% colnames(meta)){stop(paste0("signature.x: ",signature.x," is not in data,please check !"))}
  if (!group.by %in% colnames(meta)){stop(paste0("group.by: ",group.by," is not in data,please check !"))}
  data <- rds@meta.data[,c(group.by,signature)]
  colnames(data) <- c("y","x")
  font1 <- 10
  textSize <- 10
  p1 <- ggplot2::ggplot(data,aes(x=x,y=y,fill=stat(x)))+
    geom_density_ridges_gradient()+
    scale_fill_viridis_c(name = "value",option = "C") +
    labs(x=paste0("The score of ",signature),y="") +
    theme_classic()
  return(p1)
}

#' scgmt_density_plot
#'
#' @param rds      seurat object
#' @param signature
#' @param group.by
#' @param cols
#' @param alpha
#'
#' @return  ggplot2 object
#' @export
#'
#' @examples
#'
#'

scgmt_density_plot <- function(rds,signature,group.by="ident",cols, alpha=0.7){

  rds$ident <- rds@active.ident
  meta <- rds@meta.data
  if (!signature %in% colnames(meta)){stop(paste0("signature.x: ",signature.x," is not in data,please check !"))}
  if (!group.by %in% colnames(meta)){stop(paste0("group.by: ",group.by," is not in data,please check !"))}
  data <- rds@meta.data[,c(group.by,signature)]
  colnames(data) <- c("y","x")
  font1 <- 10
  textSize <- 10
  p1 <- ggplot2::ggplot(data, aes(x = x))+ geom_density(aes(fill = y,color =y), alpha=alpha) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    guides(color="none") +
    labs(x=paste0("The score of ",signature),y="Density",fill="") +
    theme_classic()
  return(p1)
}
