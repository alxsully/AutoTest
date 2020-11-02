#'@title PosLR
#'@description Finds differentially expressed genes for each cluster within a Seurat object
#'@import Seurat
#'@import tidyverse
#'@details Returns a table containing genes most specific to each cluster, including sensitivity, specificity, CI and LR for each gene
#'@param object A Seurat object
#'@param sens.min Minimum sensitivity threshold. Default is 0.50
#'@param spec.min Minimum specificity threshold. Default is 0.50
#'@param top.n.genes Number of genes returned for each cluster. Default is 10
#'@param alpha Significance level used for calculating CI. Default is 0.01
#'@param slot Slot to pull data from. Default is "data"
#'@param assay Assay to use in differential expression testing. Default is "SCT"
#'@return data.frame containing the top.n.genes most specific to each cluster
#'@examples
#' #Find top 5 DE genes for each cluster
#' PosLR(object = pbmc small, top.n.genes = 5)
#'@export
PosLR <- function(object, 
                  sens.min = 0.50,
                  spec.min = 0.50,
                  top.n.genes = 10,
                  alpha = 0.01,
                  slot = "data",
                  assay = "SCT", ...){
  # Determine dataset size
  N <- length(object$nCount_RNA)
  
  # Pull variable features
  genes.var <- VariableFeatures(object)
  
  # Create empty lists and tables for storing test results
  n <- c()
  u <- c()
  alpha.adj <- c()
  se <- ratio <- spec <- sens <- b <- a <- matrix(ncol = length(levels(object)), 
    nrow = length(genes.var))
  temp.object <- vector(mode = "list", length = length(levels(object)))

  #Progress Update
  print("Calculating sensitivity of each gene for each cluster.")

  #Function for assay data extraction
  pull.matrix <- function(x) {
    GetAssayData(object, slot = slot, assay = assay)[genes.var, WhichCells(object, ident = x)] %>%
      as.data.frame() %>%
      t()
  }
  
  #Extract assay data, stored as list
  temp.object <- lapply(levels(object), pull.matrix)
    
  # Determine and store cluster size
  n <- sapply(temp.object, nrow)
  
  #Calculate and store number of cells expressing
  a <- sapply(temp.object, FUN = function(x) {colSums(x != 0)})
  
  #Calculate and store "Sensitivity" (cells expressing/total cells)
  sens <- t(t(a)/n)
  
  #Progress Update
  print("Calculating specificity of each gene for each cluster.")
  
  #Calculate false positives
  b <- sapply(1:ncol(a), FUN = function (x) {rowSums(a[,-x])})
  
  #Calculate Specificity (# Cells not in cluster - False positives, divided by # Cells not in cluster)
  u <- N - n
  spec <- t((u-t(b))/u)
  
  #Progress Update
  print("Calculating +LR.")
  
  #Calculate and store LR+
  ratio <- sens/(1-spec)
  
  #Calculate standard error (note that sweep simply subtracts a vector
  #from each row (MARGIN =2) or column (MARGIN = 1) of a matrix) 
  se <- sqrt(sweep(1/a, MARGIN = 2, 1/n) + sweep(1/b, MARGIN = 2, 1/u))
  
  #Counts number of genes that reach the sensitivity and specificity threshold,
  #used for calculating Bonferroni-corrected alpha
  pos.results <- (spec > spec.min) & (sens > sens.min)
  count <- colSums(pos.results)
  
  # Calculate Bonferroni-corrected alpha
  alpha.adj <- alpha / count

  #Calculate and store confidence intervals
  margin.of.error <- t(qnorm(1 - alpha.adj/2) * t(se))
  lo.bound <- exp(log(ratio[pos.results]) - margin.of.error[pos.results])
  up.bound <- exp(log(ratio[pos.results]) + margin.of.error[pos.results])
  CI <- paste0("(", round(lo.bound, 2), ", ", round(up.bound, 2), ")")
  
  #Assemble results table from the following: cluster label, gene name, sensitivity, specificity, confidence interval, and LR
  test.res <- data.frame(mapply(x = levels(object), y = count, 
                                FUN = function(x,y) {rep(x,y)}, USE.NAMES = FALSE) %>%
                           unlist() %>%
                           as.numeric(),
                    apply(pos.results, MARGIN = 2,
                          FUN = function(x) {genes.var[x]}) %>% 
                      unlist(),
                    round(sens[pos.results],2),
                    round(spec[pos.results],2), 
                    CI,
                    round(ratio[pos.results],2),
                    lo.bound)
  
  #Filter for lo.bound > 1, remove lower bound column
  test.res <- test.res[test.res[,7] > 1, 1:6]

  # Format and export output table
  colnames(test.res) <- c("cluster",
                          "gene",
                          "sens",
                          "spec",
                          paste0("adj ", toString((1 - alpha) * 100), "% CI"), 
                          "LR")
  
  test.res <- test.res %>%
    dplyr::group_by(cluster) %>%
      top_n(top.n.genes) %>%
        arrange(cluster, desc(LR))

  return(test.res)
}
