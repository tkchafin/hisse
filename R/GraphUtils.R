#utility functions needed for running hisse on networks
read.snaq <- function(file){
  #NOTE: Function will only read first network in a file (e.g., from the .networks file this is the max net)
  contents <- scan(file = file, what = "", sep = ";", quiet=TRUE)
  net <- paste0(contents[1], ";")
  #move lambda values to edge labels
  net<-gsub("(#H[0-9]*?):([0-9]+?.[0-9]*?)::([0-9].[0-9]*)([\\,,\\),\\(])", "\\1_\\3:\\2\\4", net) 

  #read as phylo
  x<-read.tree(text=net, keep.multi=FALSE)

  #modified from ape::as.evonet.phylo
  pos <- grep("#", x$tip.label)
  ind <- match(pos, x$edge[, 2])
  reticulation <- x$edge[ind, , drop = FALSE]
  edge <- x$edge[-ind, , drop = FALSE]
  nTips <- as.integer(length(x$tip.label))
  
  #capture lambdas and fix node/tip labels
  gammas<-vector()
  index<-1
  for (i in 1:length(x$node.label)){
    gammas[index]<-gsub(".*?_(.*)", "\\1", x$node.label[i])
    if (gammas[index] == "") gammas[index] <- 1.0
    x$node.label[i] <- gsub("(.*?)_.*", "\\1", x$node.label[i])
    index<-index+1
  }
  for (k in pos){
    gammas[index]<-gsub(".*?_(.*)", "\\1", x$tip.label[k])
    x$tip.label[k] <- gsub("(.*?)_.*", "\\1", x$tip.label[k])
    index<-index+1
  }
  
  reticulation[, 2] <- as.integer(match(x$tip.label[pos], x$node.label) + nTips)
  for (i in sort(pos, TRUE)) {
    edge[edge > i ] <- edge[edge > i] - 1L
    reticulation[reticulation > i] <- reticulation[reticulation > i] - 1L
  }
  x$edge <- edge
  x$reticulation <- reticulation
  #
  #NOTE: Change from the stock 'evonet' object from ape -- Added 
  #a member variable storing lambda values
  #
  x$anc.prop <- gammas
  
  if (!is.null(x$edge.length)){
    x$hybrid.edge.length <- as.vector(x$edge.length[ind])
    x$edge.length <- x$edge.length[-ind]
  }
  x$tip.label <- x$tip.label[-pos]
  class(x) <- c("evonet", "phylo")
  x
}

HasCycle <- function(A){
  #this follows stackoverflow answer by user Casteels here:
  #https://stackoverflow.com/questions/16436165/detecting-cycles-in-an-adjacency-matrix/25537032#25537032
  
  #Laplacian matrix; value i,j = sum(A[i,]) if i==j; -1 if i!=j & A[i,j] == 1; else 0
  L <- -A
  for (i in 1:nrow(A)) L[i,i] <- sum(A[i,])
  
  #graph has cycles if 0.5 * trace.L >- rank.L + 1 ... don't ask me why this works.
  rank.L <- sum(A)-1 #number edges - number connected components
  trace.L <- 2*sum(A)  #2*number of edges
  if(0.5*(trace.L) >= rank.L+1){
    return(TRUE) #has cycles
  }else{
    return(FALSE) #has no cycles
  }
}

#Recursively searches a graph (in form of adjacency matrix)
#inwards towards the root from a given tip node and looks 
#for circular dependency
HasCircularDependency<- function(A){
  visited <- rep(FALSE, nrow(A))
  recStack <- rep(FALSE, nrow(A))
  for (i in 1:nrow(A)){
    if (rowSums(A)[i] != 0){ #skip terminal nodes; no point
      print(i)
      if (DepthFirstCheckCycle(A, i, visited, recStack)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

#DFS of a given sub-graph, looking for circular dependencies
DepthFirstCheckCycle <- function(A, i, visited, recStack){
  if (visited[i] == FALSE){
    visited[i] <- TRUE
    recStack[i] <- TRUE
    
    #recursion for all adjacent nodes
    for (j in which(A[i,]==1)){
      if (!visited[j]){
        if (DepthFirstCheckCycle(A, j, visited, recStack)){
          return(TRUE)
        }
      }else if (recStack[j] == TRUE){
        return(TRUE)
      }
    }
  }
  recStack[i] <- FALSE
  return(FALSE)
}

#Takes an ape::evonet object and builds an adjacency matrix
#technically would also work on ape::phylo but not sure if 
#there would ever be a reason to use it
GetAdjacencyMatrix <- function(net){
  n <- net$Nnode+length(net$tip.label)
  adj <- matrix(0, n, n)
  
  #for each edge, e[1] depends on e[2]
  for (e in 1:nrow(net$edge)){
    adj[net$edge[e,1], net$edge[e,2]] <- 1
  }
  
  #for reticulate edges, r[2] depends on r[1]
  if (!is.null(net$reticulation)){
    for (r in 1:nrow(net$reticulation)){
      adj[net$reticulation[r,2], net$reticulation[r,1]] <- 1
    }
  }
  return(adj)
}

#function ranks nodes by neighborhood size, returning
#a list of vectors where names(list) = neighborhood sizes (e.g., 
#1 = terminal nodes and max(names(list)) = root
#NOTE that 'empty' neighborhood size groups are removed
RankByNeighborhoodSize <- function(A){
  sizes <- NeighborhoodSizes(A)
  ranks <- rep(list(list()), (max(sizes)+1))
  for (i in 1:(max(sizes)+1)){
    indices<-which(sizes==(i-1))
    if (length(indices)>0){
      #print(indices)
      ranks[i]<-list(indices)
    }
  }
  names(ranks) <- 1:length(ranks)
  ranks <- ranks[lengths(ranks)!=0]
  return(ranks)
}

#function calculates total descendents for nodes in a 
#directed graph given an adjacency matrix
NeighborhoodSizes <- function(A){
  sizes <- rep(0, nrow(A))
  #for each internal node, recurse through descendants 
  #and increment value in the 'depth' table
  for (i in which(rowSums(A)!=0)){
    #print(paste0("Node: ",i))
    sizes[i] <- NeighborhoodSizeFocalNode(A, i, sizes)
  }
  return(sizes)
}

NeighborhoodSizeFocalNode <- function(A, i, sizes){
  total<-0
  #if node has already been visited, return its size
  if (sizes[i] != 0 || rowSums(A)[i]==0){
    return(sizes[i])
  }
  #else, recurse over descendants
  for (j in which(A[i,]==1)){
    #print(paste0("  - Descendnt: ",j))
    total <- total+1
    total <- total+(NeighborhoodSizeFocalNode(A, j, sizes))
  }
  sizes[i] <- total
  return(total)
}