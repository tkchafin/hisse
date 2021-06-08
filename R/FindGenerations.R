## Common function used by multiple methods.

#Note that James Boyko came up with the FindGenerations code, so the insanity of this Rumsfeldian speak is all him....
FindGenerations <- function(phy){
    generation <- list()
    known <- 1:Ntip(phy)
    unknown <- phy$edge[,1]
    needed <- phy$edge[,2]
    root <- min(unknown)
    i <- 1
    repeat{
        knowable <- unknown[needed %in% known]
        knowable <- knowable[duplicated(knowable)]
        generation[[i]] <-  knowable
        
        known <- c(known, knowable)
        needed <- needed[!unknown %in% knowable]
        unknown <- unknown[!unknown %in% knowable]
        i <- i + 1
        if (any(root == knowable)) break
    }
    res <- generation
    return(res)
}

FindGenerations.net <- function(net){
  #build adjacency matrix 
  A <- GetAdjacencyMatrix(net)
  
  #check for cycles; if TRUE check if they involve circular dependency
  print("checking for cycles")
  if (HasCycle(A)){
    if (HasCircularDependency(A)){
      #NOTE: This function will only return index of the first circular dependency it finds.
      stop(paste0("Network has circular dependency; remove it before proceeding."))
    }
  }
  
  #rank nodes by number of dependencies
  ranks<-RankByNeighborhoodSize(A)
  names(ranks) <- NULL
  ranks <- ranks[2:length(ranks)]

  #TO-DO  
  #Move nodes downwards to balance gen sizes
  #(this is so we can parallelize later if we want)
  #ranks<-OptimizeRanks(ranks)
}

#deprecated, version based on OG FindGenerations
#problem with this version is that it is infinite if 
#the root is never knowable -- need to check for 
#circular dependencies
# FindGenerations.net <- function(phy){
#   generation <- list()
#   known <- 1:Ntip(phy)
#   unknown <- c(phy$edge[,1], phy$reticulation[,2])
#   needed <- c(phy$edge[,2], phy$reticulation[,1])
#   root <- min(unknown)
#   i <- 1
#   repeat{
#     #
#     #Only real change here is that I added the nodes
#     #from the reticulation table (phy$reticulation) to knowable
#     # -- this allows originations on terminal branches 
#     #to be included, otherwise they get missed in the first
#     #iteration here
#     knowable <- unknown[needed %in% known]
#     knowable <- c(knowable, as.vector(phy$reticulation[,1]))
#     knowable <- knowable[duplicated(knowable)]
#     generation[[i]] <-  knowable
# 
#     known <- c(known, knowable)
#     needed <- needed[!unknown %in% knowable]
#     unknown <- unknown[!unknown %in% knowable]
#     i <- i + 1
#     if (any(root == knowable)) break
#   }
#   res <- generation
#   return(res)
# }
