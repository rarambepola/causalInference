##function that orients j - k to j -> k whenever exists i -> j st i and k
##are not adjacent (otherwise this would be a v-structure)
.rule1 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(dirEdge(adjmat, i, j)){
        for(k in (1:nvar)[-i]){
          if(undirEdge(adjmat, j, k) && !anyEdge(adjmat, i, k)){
            adjmat[k,j] <- 0
            changed <- TRUE
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}

##function that orients i-j to i->j whenever there is i->k->j (otherwise there
##would be a cycle)
.rule2 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(k in 1:nvar){
      if(dirEdge(adjmat, i, k)){
        for(j in (1:nvar)[-i]){
          if(dirEdge(adjmat, k, j) && undirEdge(adjmat, i, j)){
            adjmat[j,i] <- 0
            changed <- TRUE
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}

##function that orients i-j to i->j whenever there are chains i-k->j and i-l->j
##such that k, l are not adjacent (otherwise a new v-structure or cycle would be formed)
.rule3 <- function(adjmat){
  changed <- FALSE
  nvar <- dim(adjmat)[1]
  for(i in 1:nvar){
    for(j in 1:nvar){
      if(undirEdge(adjmat, i, j)){
        for(k in (1:nvar)[-c(i, j)]){
          if(undirEdge(adjmat, i, k) && dirEdge(adjmat, k, j)){
            for(l in (1:nvar)[-c(i, j, k)]){
              if(undirEdge(adjmat, i, l) && dirEdge(adjmat, l,l)){
                adjmat[j, i] <- 0
                changed <- TRUE
              }
            }
          }
        }
      }
    }
  }
  return(list(adjmat, changed))
}


##function to get edges from (undirected) adjacency matrix
getEdges <- function(adjMat){
  n <- dim(adjMat)[1]
  edges <- list()
  for(i in 1:n){
    for(j in 1:n){
      if(adjMat[i, j] == 1){
        edges[[length(edges) + 1]] <- sort(c(i, j))
      }
    }
  }
  return(unique(edges))
}

##function to get set Adj(G, i)\j
getAdj <- function(adjMat, i, j){
  n <- dim(adjMat)[1]
  adjs <- c()
  for(k in (1:n)[-j]){
    if(adjMat[k, i] == 1){
      adjs <- c(adjs, k)
    }
  }
  return(adjs)
}

##function to make complete adjacency matrix
compltMat <- function(n){
  return(matrix(1, nrow=n, ncol=n) - diag(n))
}

##function that takes a set and a number and returns all
##subsets of that size
getSubsets <- function(set, nSubset){
  sslist <- list()
  if(nSubset == 1){
    for(i in 1:length(set)){
      sslist[[i]] <- c(set[i])
    }
  }else{
    for(i in 1:(length(set) - nSubset + 1)){
      sslist2 <- getSubsets(set[-(1:i)], nSubset - 1)
      for(j in 1:length(sslist2)){
        sslist[[length(sslist) + 1]] <- c(set[i], sslist2[[j]])
      }
    }
  }
  return(sslist)
}


##function that takes an adjacency matrix and returns
##unshielded triples
getUnTrpls <- function(adjmat){
  triples <- list()
  nvars <- dim(adjmat)[1]
  for(i in 1:nvars){
    for(j in 1:nvars){
      if(adjmat[i, j] == 1){
        for(k in (1:nvars)[-i]){
          if(adjmat[j, k] == 1){
            if(adjmat[i, k] == 0){
              if(i < k){
                triples[[length(triples) + 1]] <- c(i, j, k)
              }else{
                triples[[length(triples) + 1]] <- c(k, j, i)
              }
            }
          }
        }
      }
    }
  }
  triples <- unique(triples)
  return(triples)
}

##functions that returns boolean if:
##there is an edge i -> j
dirEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 && adjmat[j, i] == 0)
}
##there is an edge i <->j
undirEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 && adjmat[j, i] == 1)
}
##there is an edge i - j
anyEdge <- function(adjmat, i, j){
  return(adjmat[i, j] == 1 || adjmat[j, i] == 1)
}

