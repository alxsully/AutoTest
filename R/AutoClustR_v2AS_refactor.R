#AutoClustR
#'@title AutoClustR
#'@description Determines the optimal number of principal components and optimal values for the resolution and K.param arguments used for clustering
#'@import tidyverse
#'@import Seurat
#'@import segmented
#'@import cowplot
#'@import clusterCrit

#'@param object The Seurat object to be tested
#'@param method Specify method used for cluster ranking. "CH" for Calisnki-Harabasz,
#' "DB" for Davies-Bouldin, "Sil" for silhouette. Default is "CH"
#'@param standardize Determines whether or not to standardize principal components. Default is TRUE
#'@param k.param.space Numeric vector of length three, specifying minimum, maximum, and step size. Default is c(10, 240, 10)
#'@param resolution.space Numeric vector of length three, specifyin gminimum, maximum, and step size. Default is c(0.1, 2.4, 0.1)
#'@param n.iterations Number of restarts for hill-climbing optimization.
#'Lower values will give lower runtimes but may not find optimal parameter values. Default is 24
#'@param pc.use Number of principal components used for segmented linear regression testing. Default is 200
#'@param pc.estimate Used by segemented function to estimate segmented linear regression breakpoint. Default is 10.
#'@return A list containing a seurat object clustered using optimal parameters, a double containing the number of PC's used,
#'two matrices containing the number of clusters and the ranking for tested parameter combinations,
#'an elbow plot showing the segmented linear regression, and a heatmap showing ranknigs for tested parameter combinations
#'@examples
#' #Returns list with clustered seurat object and descriptive data
#' clustering_data <- AutoClustR(seurat.object, method = "CH", k.param.space = c(10, 120, 10))
#'@export
AutoClustR <- function(object,
                       method = "CH",
                       standardize = TRUE,
                       k.param.space = c(10, 240, 10),
                       resolution.space = c(0.1, 2.4, 0.1),
                       n.iterations = 24,
                       pc.use = 200,
                       pc.estimate = 10, ...){
  # Step 1: Format inputs and screen for errors
  CheckSpace(k.param.space, resolution.space)
  method <- CheckMethod(method)

  # Load principal components and singular values
  object <- CheckPCs(object, pc.use = pc.use)

  # Step 2: Use segmented regression to determine psi (point of intersection for linear regression.)
  model <- ChooseDimensionality(object, pc.use, pc.estimate)
  npc <- floor(model$psi[[2]])

  #step 3: Extract relevant PCs for all cells for use in subsequent analyses
  embeddings <- LoadPCA(object, npc = npc, standardize = standardize)

  #Step 4: Construct the parameter space to be searched
  param.space <- ConstructParamSpace(k.param.space, resolution.space)

  #Step 5: Use random restart hill climbing to find local maxima/minima of internal clustering validation index
  hill.results <- HillClimb(object,
                            param.space = param.space,
                            n.iterations = n.iterations,
                            embeddings = embeddings,
                            npc = npc,
                            method = method)
  #Step 6: Construct relevant plots, package results as list for export
  results <- ConstructResults(object,
                              best.params = hill.results[[1]],
                              indices = hill.results[[2]],
                              n.clusters = hill.results[[3]],
                              pc.use = pc.use,
                              npc = npc,
                              model = model,
                              method = method,
                              standardize = standardize)
  return(results)
}

is.whole <- function(x){
  if((x %% 1) == 0){TRUE}
  else{FALSE}
}

is.divisible <- function(space){
  round(((space[2] - space[1]) / space[3]), digits = 5) %>%
    is.whole()
}

# Test if step sizes divide evenly into paramter space
CheckSpace <- function(k.param.space, resolution.space) {
  if(!is.whole(k.param.space[3])){
    stop("k.param step size must be a whole number")
  } else if(!is.divisible(k.param.space)){
    stop("K.param range is not divisible by step size")
  } else if(!is.divisible(resolution.space)){
    stop("Resolution range is not divisible by step size")
  }
}

# Format 'method' for clusterCrit
CheckMethod <- function(method){
  if(method == "CH"){
    method <- "Calinski_Harabasz"
  } else if(method == "DB"){
    method <- "Davies_Bouldin"
  } else if(method == "Sil"){
    method <- "Silhouette"
  } else{
    stop("Index not supported.")
  }
}

#Write additional PCs (if necessary)
CheckPCs <- function(object, pc.use){
  if(length(object[["pca"]]@stdev) < pc.use){
    print("Calculating additional principal components")
    object <- RunPCA(object, npcs = pc.use, verbose = FALSE)
  }
  return(object)
}

#Use segmented regression to determine psi (point of intersection for linear regression.)
ChooseDimensionality <- function(object, pc.use, pc.estimate){
  # Perform segmented regression
  pcs <- 1:pc.use
  std_dev <- object[["pca"]]@stdev[pcs]
  pc.data <- data.frame(pcs, std_dev)
  print("Choosing dimensionality.")
  model <- segmented(lm(formula = std_dev ~ pcs, data = pc.data),
                     seg.Z = ~pcs,
                     psi = pc.estimate,
                     control = seg.control(seed = 1))
  return(model)
}

# Load PCA results
LoadPCA <- function(object, npc, standardize){
    all.embeddings <- Embeddings(object)[, 1:npc]
    std_dev <- object[["pca"]]@stdev[1:npc]
    # Standardize principal components if applicable
    # Theres probably a better way to do this using apply vs a for loop
    if(standardize == TRUE){
      for(i in 1:npc){
        all.embeddings[, i] <- (all.embeddings[, i] - mean(all.embeddings[, i])) / std_dev[i]
        }
    }
    return(all.embeddings)
}

ConstructParamSpace <- function(k.param.space, resolution.space){
  param.space <- list(seq(k.param.space[1], k.param.space[2], k.param.space[3]),
                      seq(resolution.space[1], resolution.space[2], resolution.space[3]))
  names(param.space) <- c("k.param", "resolution")
  return(param.space)
}

#Excessively Long function that calculates index values for parameter pairs
HillClimb <- function(object, param.space, n.iterations, embeddings, npc, method){

  #Use Latin square sampling to randomize starting parameters
  latin.square <- LatinSquare(param.space, n.iterations)

  # Create empty tables for storing AutoClustR results
  #This would be better done in such a way that indices and n.clusters are part of the same object,
  #especially because they have the same dim names
  indices <- n.clusters <- matrix(0,
                                  nrow = length(param.space$k.param),
                                  ncol = length(param.space$resolution),
                                  dimnames = list(K.param = param.space$k.param,
                                                  Resolution = param.space$resolution))

  # Create progress bar
  print("Choosing clustering parameters.")
  prog.bar <- txtProgressBar(min = 0, max = n.iterations, char = "+", style = 3)

  # (Re)start hill climbing
  for(m in 1:n.iterations){
    #Set start point for iteration m
    i <- latin.square$i[m]
    j <- latin.square$j[m]

    #Perform Cluster analysis on start point
    if (index.is.valid(i, j, indices)){
      results          <- AnalyzeClusters(i, j, object,
                                          npc = npc,
                                          param.space = param.space,
                                          method = method,
                                          embeddings = embeddings)
      indices[i, j]    <- results[1]
      n.clusters[i, j] <- results[2]
    }

    #Search adjacent parameter pairs until local maximum is found
    while(adjacent.are.zero(i, j, indices)){
      # Perform cluster analysis for parameter pairs adjacent to start point
      # ie [i - 1, j], [i + 1, j], [i, j - 1], [i, j + 1]
      if (index.is.valid(i - 1, j, indices)) {
        results              <- AnalyzeClusters(i - 1, j, object,
                                                npc = npc,
                                                param.space = param.space,
                                                method = method,
                                                embeddings = embeddings)
        indices[i - 1, j]    <- results[1]
        n.clusters[i - 1, j] <- results[2]
      }
      if (index.is.valid(i + 1, j, indices)) {
        results              <- AnalyzeClusters(i + 1, j, object,
                                                npc = npc,
                                                param.space = param.space,
                                                method = method,
                                                embeddings = embeddings)
        indices[i + 1, j]    <- results[1]
        n.clusters[i + 1, j] <- results[2]
      }
      if (index.is.valid(i, j - 1, indices)) {
        results              <- AnalyzeClusters(i, j - 1, object,
                                                npc = npc,
                                                param.space = param.space,
                                                method = method,
                                                embeddings = embeddings)
        indices[i, j - 1]    <- results[1]
        n.clusters[i, j - 1] <- results[2]
      }
      if (index.is.valid(i, j + 1, indices)) {
        results              <- AnalyzeClusters(i, j + 1, object,
                                                npc = npc,
                                                param.space = param.space,
                                                method = method,
                                                embeddings = embeddings)
        indices[i, j + 1]    <- results[1]
        n.clusters[i, j + 1] <- results[2]
      }

      # Determine next step for hill climbing
      next.step <- FindNextStep(i, j, indices, method)
      i <- next.step$i
      j <- next.step$j
    }

    # Update progress bar
    setTxtProgressBar(prog.bar, m)
  }

  # Input 'NA' for untested solutions
  n.clusters[n.clusters == 0] <- NA
  indices[indices == 0] <- NA

  # Determine 'best' parameters
  if(method == "Davies_Bouldin"){
    best.index <- which(indices == min(indices, na.rm = TRUE), arr.ind = TRUE)
  } else{
    best.index <- which(indices == max(indices, na.rm = TRUE), arr.ind = TRUE)
  }

  best.params <- list(param.space$k.param[[best.index[1, 1]]],
                      param.space$resolution[[best.index[1, 2]]])
  names(best.params) <- c("k", "res")
  return(list(best.params, indices, n.clusters))
}

ConstructResults <- function(object, best.params, indices, n.clusters, pc.use, npc, model, method, standardize){
  # Construct SNN graph for ("K.param", "Resolution") = (best.K.param, best.Resolution)
  object <- FindNeighbors(object,
                          k.param = best.params$k,
                          compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                          nn.eps = 0.0, verbose = FALSE, reduction = "pca", dims = 1:npc) %>%

    # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (best.K.param, best.Resolution)
    FindClusters(modularity.fxn = 1,
                 resolution = best.params$res,
                 algorithm = 1,
                 group.singletons = TRUE,
                 verbose = FALSE)

  # Format and export segmented regression plot
  pcs <- 1:pc.use
  m <- slope(model)$pcs
  b <- intercept(model)$pcs
  psi <- model$psi[[2]]
  std.dev <- object[["pca"]]@stdev[pcs]

  npc.plot <- ggplot(as.data.frame(cbind(pcs, std.dev)), aes(x = pcs, y = std.dev)) +
    geom_point(size = 2.5) +
    geom_segment(aes(x = 1, y = (m[1, 1] * 1) + b[1, 1], xend = psi,
                     yend = (m[1, 1] * psi) + b[1, 1]), size = 1) +
    geom_segment(aes(x = psi, y = (m[2, 1] * psi) + b[2, 1], xend = pc.use,
                     yend = (m[2, 1] * pc.use) + b[2, 1]), size = 1) +
    geom_vline(xintercept = npc, color = "red3", size = 1, linetype = "longdash") +
    geom_text(aes(x = npc + (pc.use / 20), label = paste0("npc = ", npc),
                  y = max(std.dev, na.rm = TRUE) / 2), color = "red3", size = 5) +
    theme_cowplot(font_size = 16) +
    xlab("PC") +
    ylab("Standard Deviation")

  # Format and export AutoClustR heatmap
  clust.tab <- as.data.frame.table(n.clusters)
  ind.tab <- as.data.frame.table(indices)
  index <- paste0(gsub("[^::A-Z::]","", method), "I")
  index <- if(standardize == TRUE){paste0(index, "_Std")}

  for(i in 1:length(ind.tab[, 1])){
    if(ind.tab[i, 1] == best.params$k && ind.tab[i, 2] == best.params$res){
      clust.tab[i, 4] <- ind.tab[i, 4] <- "X"
    }else{
      clust.tab[i, 4] <- ind.tab[i, 4] <- ""
    }
  }

  colnames(clust.tab) <- colnames(ind.tab) <- c("K.param", "Resolution", index, "best.index")

  ind.plot <- ggplot(ind.tab, aes_string(x = "K.param", y = "Resolution", fill = index)) +
    geom_tile() +
    scale_fill_gradient(low = "gold1", high = "darkred") +
    geom_text(aes_string(label = "best.index"), size = 6) +
    theme(axis.title = element_text(size = 18),
                                                                                                                                                                                                               axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 15),
                                                                                                                                                                                                               legend.text = element_text(size = 15), legend.title = element_text(size = 18)) + labs(x = "\nK.param", y = "Resolution\n", fill = index)

  # Return outputs
  return(list(object, npc, n.clusters, indices, npc.plot, ind.plot))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Ancillary Functions for Hill Climb

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LatinSquare <- function(param.space, n.iterations){
  set.seed(1)
  i0 <- c()
  j0 <- c()

  square.size <- c(length(param.space$k.param), length(param.space$resolution)) / n.iterations
  sample <- list(sample(1:n.iterations, n.iterations, replace = FALSE),
                 sample(1:n.iterations, n.iterations, replace = FALSE))

  if(sum(square.size %% 1) != 0) {
    stop("Error: Number of steps in parameter space must be evenly divisible by n.iterations")
    #If it's not evenly divisible, you get errors trying to pass non-integer values
    #As indices to tables. However, there's probably an easy way to fix this.
  }

  i0 <- square.size[1] * (sample[[1]] - 1) + sample(1:square.size[1], n.iterations, replace = TRUE)
  j0 <- square.size[2] * (sample[[2]] - 1) + sample(1:square.size[2], n.iterations, replace = TRUE)
  latin.square <- list(i0,j0)
  names(latin.square) <- c("i", "j")
  return(latin.square)
}

index.is.valid <- function(i,j, indices){
  if (i > 0 &
      j > 0 &
      i  <= nrow(indices) &
      j  <= ncol(indices)){
    if(!is.na(indices[i,j])) {
      if(indices[i, j] == 0){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

ConstructNeighbors <- function(i, j, indices){
  neighbors <- matrix(c(
    if((i - 1) > 0){indices[i - 1, j]} else{NA},
    if((j - 1) > 0){indices[i, j - 1]} else{NA},
    if((i + 1) <= nrow(indices)){indices[i + 1, j]} else{NA},
    if((j + 1) <= ncol(indices)){indices[i, j + 1]} else{NA}),
    nrow = 2, ncol = 2)
  return(neighbors)
}

adjacent.are.zero <- function(i, j, indices){
  neighbors <- ConstructNeighbors(i, j, indices)
  are.zero <- suppressWarnings(min(neighbors, na.rm = TRUE) == 0)
  return(are.zero)
}

#Returns index value and number of clusters for given parameter pair (i,j)
AnalyzeClusters <- function(i, j, object, npc, param.space, method, embeddings){
    # Construct SNN graph for ("K.param", "Resolution") = (i, j)
    # These need to be stored. As it is, the nn graph is being calculated multiple times for identical k.params
    temp.object <- FindNeighbors(object,
                                 k.param = param.space$k.param[i],
                                 compute.SNN = TRUE,
                                 prune.SNN = 1/15,
                                 nn.method = "rann",
                                 nn.eps = 0.0,
                                 verbose = FALSE,
                                 reduction = "pca",
                                 dims = 1:npc) %>%

      # Cluster cells by maximizing modularity for ("K.param", "Resolution") = (i, j)
      FindClusters(modularity.fxn = 1,
                   resolution = param.space$resolution[j],
                   algorithm = 1,
                   group.singletons = TRUE,
                   verbose = FALSE)

    # Determine and store cluster number for ("K.param", "Resolution") = (i, j)
    clusters <- as.integer(Idents(temp.object))
    K <- length(levels(temp.object))

    # Calculate and store internal clustering validation index for ("K.param", "Resolution") = (i, j)
    if(K == 1){
      ind <- NA
    } else{
      ind <- clusterCrit::intCriteria(embeddings, part = clusters, crit = method)[[1]]
    }
  return(c(ind, K))
}

FindNextStep <- function(i, j, indices, method){
  if(is.na(indices[i, j])){
  next.step <- list(i, j)
  names(next.step) <- c("i", "j")
  return(next.step)
  } else {
    if(method == "Davies_Bouldin"){
      next.step <- FindBestNeighbor(i, j, indices = -indices)
      return(next.step)
    } else{
      next.step <- FindBestNeighbor(i, j, indices)
      return(next.step)
    }
  }
}

FindBestNeighbor <- function(i, j, indices) {
  neighbors <- ConstructNeighbors(i, j, indices)
  if(indices[i, j] >= max(neighbors, na.rm = TRUE) || all(is.na(neighbors))){
    i <- i
    j <- j
  } else{
    best.neighbor <- which(neighbors == max(neighbors, na.rm = TRUE),
                           arr.ind = TRUE)
    if(best.neighbor[1, 1] == 1){
      i <- i - 3 + (2 * best.neighbor[1, 2])
      j <- j

    } else{
      i <- i
      j <- j - 3 + (2 * best.neighbor[1, 2])
    }
  }
  next.step <- list(i,j)
  names(next.step) <- c("i","j")
  return(next.step)
}

