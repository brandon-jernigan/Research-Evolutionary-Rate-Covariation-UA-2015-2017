#!/usr/bin/env Rscript
#Script pulled from a larger script provided by Dr. Nathan Clark

library(methods)

if(!require("ape")){
  install.packages("ape", repos="http://cran.r-project.org")
}
if (!require("compiler")){
  install.packages("compiler", repos="http://cran.r-project.org")
}
if(!require("phytools")){
  install.packages("phytools", repos="http://cran.r-project.org")
}

getRoot = function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]

getChildren=function(tree, nodeN){
  tree$edge[tree$edge[,1]==nodeN,2]
}

getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }
  
}

treeTraverse=function(tree, node=NULL){
  if(is.null(node)){
    rt=getRoot(tree)
    ic=getChildren(tree,rt)
    return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
    
  }
  else{
    if (node<=length(tree$tip)){
      return(tree$tip[node])
    }
    else{
      ic=getChildren(tree,node)
      return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))
      
    }
  }
}

#puts tree it tip sorted order
CanonicalForm=function(tree){
  par(mfrow=c(1,2))
  oo=order(tree$tip.label)
  tree$tip.label=tree$tip.label[oo]
  ii=match(1:length(oo), tree$edge[,2])
  tree$edge[ii,2]=order(oo)
  rotateConstr(tree, sort(tree$tip.label))
  
}

#computes all paths from the tips to ancestors
allPaths=function(tree){
  dd=dist.nodes(tree)
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      allD=c(allD, dd[i, ia])
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))
    }
  }
  return(list(dist=allD, nodeId=nn))
}


#mathces node if descendands in tr1 are a subset of descendands in tr2
matchNodesInject=function (tr1, tr2){
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  for (i in 1:length(desc.tr1)) {
    Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
    for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                          desc.tr2[[j]]))
      Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
  }
  
  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  Nodes
}
matchNodesInject_c=cmpfun(matchNodesInject)


matchAllNodes=function(tree1, tree2){
  map=matchNodesInject_c(tree1,tree2)
  map=map[order(map[,1]),]
  map
}
matchAllNodes_c=cmpfun(matchAllNodes)


#computes paths relative to a master tree, returns a vector that is the same length as the number of
#paths in the master tree though with missing values
allPathMasterRelative=function(tree, masterTree, masterTreePaths=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPaths(masterTree)
  }
  nnMaster=masterTreePaths$nodeId[,1]*1000+masterTreePaths$nodeId[,2]
  
  treePaths=allPaths(tree)
  
  map=matchAllNodes_c(tree,masterTree)
  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]
  
  nnTree=treePaths$nodeId[,1]*1000+treePaths$nodeId[,2]
  ii=match(nnTree, nnMaster)
  #show(nnTree)
  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  vals
  
}

#a hack, turns a 2-tuble into a single value for indexing
namePaths=function(nodeMat, invert=F, mult=1000){
  if(invert){
    nodeMat=nodeMat[,c(2,1)]
  }
  return(nodeMat[,1]*mult+nodeMat[,2])
}

readTrees=function(file, max.read=NA, rearrange=F, computePaths=F){
  tmp=scan(file, sep="\t", what="character")
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  treenames=character()
  maxsp=0; # maximum number of species
  #show(length(trees))
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=unroot(read.tree(text=tmp[i]))
      #check if it has more species
      if(length(trees[[i/2]]$tip.label)>maxsp){
        maxsp=length(trees[[i/2]]$tip.label)
        allnames=trees[[i/2]]$tip.label
      }
    }
    
  }
  names(trees)=treenames
  treeObj=vector(mode = "list")
  treeObj$trees=trees
  treeObj$numTrees=length(trees)
  treeObj$maxSp=maxsp
  
  message(paste("max is ", maxsp))
  
  report=matrix(nrow=treeObj$numTrees, ncol=maxsp)
  colnames(report)=allnames
  
  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)
    
  }
  treeObj$report=report
  
  
  
  ii=which(rowSums(report)==maxsp)
  
  #Create a master tree with no edge lengths
  master=trees[[ii[1]]]
  # write.tree(master, file = "trees/master_tree_13")
  master$edge.length[]=1
  treeObj$masterTree=master
  
  
  if(rearrange){
    show(treeObj$masterTree$tip)
    treeObj$masterTree=rotateConstr(treeObj$masterTree, sort(treeObj$masterTree$tip.label))
    #this gets the abolute alphabetically constrained order when all branches
    #are present
    tiporder=treeTraverse(treeObj$masterTree)
    show(tiporder)
    #treeObj$masterTree=CanonicalForm(treeObj$masterTree)
    show(treeObj$masterTree$tip)
    for ( i in 1:treeObj$numTrees){
      
      treeObj$trees[[i]]=rotateConstr(treeObj$trees[[i]], tiporder)
      
    }
    
  }
  if (computePaths){
    ap=allPaths(master)
    #  show(length(ap))
    paths=matrix(nrow=treeObj$numTrees, ncol=length(ap$dist))
    for( i in 1:treeObj$numTrees){
      paths[i,]=allPathMasterRelative(treeObj$trees[[i]], master, ap)
    }
    colnames(paths)=namePaths(ap$nodeId)
    treeObj$paths=paths
  }
  
  treeObj
}

#drop.tip wrapper
pruneTree=function(tree, tip.names){
  keep=intersect(tree$tip.label, tip.names)
  torm=setdiff(tree$tip.label, keep)
  tree=drop.tip(tree, torm)
  tree
}

rescaleTree=function(tree){
  tree$edge.length=tree$edge.length/sqrt(sum(tree$edge.length^2))
  tree
}


matchAllNodes=function(tree1, tree2){
  map=matchNodesInject_c(tree1,tree2)
  map=map[order(map[,1]),]
  map
}

matchAllNodes_c=cmpfun(matchAllNodes)

#Maps edges to paths in a master tree
edgeIndexRelativeMaster=function(tree, masterTree){
  map=matchAllNodes(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}

scaleDist=function(x){
  x/sqrt(sum(x^2))
}
scaleDist_c=compiler::cmpfun(scaleDist)

scaleMat=function(mat){t(apply(mat,1,scaleDist_c))}

scaleMat_c=cmpfun(scaleMat)

projection <- function(protein, rna=rvector, method=c("RNA","AVE","PCA"), returnNV=F)
{
  # *************************************************************************
  # FILE: projection.R
  # AUTHOR: Tetsuya Sato <sato@kuicr.kyoto-u.ac.jp>
  # COPYRIGHT (C) 2005
  # CREATE DATE: 3/28/2005
  # MEHOD: Projection operator method
  # VERSION: 0.01
  # UPDATE DATE: 3/28/2003
  # DESCRIPTION: Given a data matrix consisting of phylogenetic vectors,
  #              compute a data matrix by removing the phylogenetic
  #              relationship among organisms applying a projection operator.
  # INPUT(arguments):
  #   1st argument: data matrix consisting of phylogenetic vectors
  #   2nd argument: phylogenetic vector based on rRNA 16S
  #   3rd argument: type of the projection
  #       (if input RNA, use RNA vector)
  #       (if input AVE, use the mean of the data matrix)
  #       (if input PCA, use the PC 1 score in the PCA analysis)
  # OUTPUT:
  #   a data matrix as a result of removing the phylogenetic relationship
  #   among organisms from the original data matrix
  # *************************************************************************
  
  
  ###projection(average vector)###
  if (match.arg(method) == "AVE") {
    n <- ncol(protein)
    m <- nrow(protein)
    protein.sc <- matrix(0,m,n)
    for (i in 1:n) {
      unit <- protein[,i] / sqrt(sum(protein[,i]^2,na.rm=TRUE))
      protein.sc[,i] <- unit
    }
    mean.clean = function(x){mm=mean(x,na.rm=TRUE); mm }
    #  average <- apply(protein.sc, 1, mean.clean)
    average=rowMeans(protein.sc, na.rm=T)
    avector.sc <- average/sqrt(sum(average^2))
    # write.table( t(avector.sc), "avg/avg_13.tsv", sep="\t", row.names=F, col.names=T , quote=F)
    #    result <- average/sqrt(sum(average^2))
    normv=avector.sc
    Id=diag(nrow(protein))
    #result=(Id-avector.sc %*% t(avector.sc))%*%protein
    result <- as.matrix(protein) - avector.sc %*% t(avector.sc) %*% as.matrix(protein)
  }
}


correlateTreesAll=function(treeObj,  usePaths=T, maxDo=NULL, species.threshold=10, species.list=NULL){
  cnt = 0
  
  maxT=treeObj$numTrees
  if (is.null(maxDo)){
    maxDo=maxT*(maxT-1)
  }
  corout=matrix(nrow=maxT, ncol=maxT)
  maxn=treeObj$report[, species.list]%*%t(treeObj$report[, species.list])
  colnames(corout)=rownames(corout)=names(treeObj$trees)
  
  done=0
  todo=length(maxn[upper.tri(maxn)]>= species.threshold)
  corout[maxn<species.threshold]=0
  diag(corout)=1
  corout[lower.tri(corout)]=0
  for (i in 1:(maxT-1)){
    for(j in (i+1):maxT){
      #show(c(i,j))
      if (is.na(corout[i,j])){
        tree1=treeObj$trees[[i]]
        tree2=treeObj$trees[[j]]
        if(!is.null(species.list)){
          tree1=pruneTree(tree1, species.list)
        }
        both=intersect(tree1$tip.label, tree2$tip.label)
        tree1=unroot(pruneTree(tree1, both))
        tree2=unroot(pruneTree(tree2, both))
        allreport=treeObj$report[,both]
        ss=rowSums(allreport)
        iiboth=which(ss==length(both))
        torm=setdiff(treeObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        if(length(both)<species.threshold){
          next
        }
        
        if(! usePaths){
          torm=setdiff(treeObj$masterTree$tip.label, both)
          allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
          for ( k in 1:length(iiboth)){
            tmptree=rescaleTree(unroot(drop.tip(treeObj$trees[[iiboth[k]]], torm)))
            allbranch[k, ]=tmptree$edge.length
          }
        }
        else{
          ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
          ii= match(namePaths(ee,T), colnames(treeObj$paths))
          allbranch=treeObj$paths[iiboth,ii]
          allbranch=scaleMat_c(allbranch)
        }
        nb=length(both)
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
        # 
        # cnt = cnt + 1
        # print(paste0("writing proj.tsv ", cnt))
        # write.table( t(proj), paste0("proj/proj" , cnt, ".tsv"), sep="\t", row.names=F, col.names=T , quote=F)		# writes tab-delimited table for other applications
        # 
        #i1=match(i, iiboth)
        #j1=match(j,iiboth)
        #corout[i,j]=cor(proj[i1, ], proj[j1,])
        #  tmpcor=cor(t(proj))
        ai=which(maxn[iiboth, iiboth]==nb, arr.ind = T)
        for (m in 1:nrow(ai)){
          k=sort(ai[m,])[1]
          l=sort(ai[m,])[2]
          
          tmpcor=cor(proj[k,], proj[l,])
          if (is.na(tmpcor)){
            tmpcor=100
          }
          corout[iiboth[k], iiboth[l]]=tmpcor
        }
        done=done+nrow(ai)
        message(paste("Done with",done, "out of", todo))
        message(paste(sum(is.na(corout))," left"), appendLF = T)
        if(done>=maxDo){
          message("Finished. Exiting because done>=maxDo")
          diag(corout)=1
          corout[maxn<species.threshold] = 100
          corout[lower.tri(corout)]=0
          corout[corout == 100] = NA
          return(corout)
        }
      }
    }
  }
  message("Finished. Exiting due to completed traversal of matrix")
  diag(corout)=1
  corout[maxn<species.threshold] = 100
  corout[lower.tri(corout)]=0
  corout[corout == 100] = NA
  corout
}

removeEmpties=function(table){
  #create vector of column/row indices to keep
  keep = apply(erc, 1, sum) != "1" | apply(erc, 2, sum) != "1"  	#the strange "1" is due to an unexplained phenomenon of the first column being summed to a nonnumeric 1. ???
  table[keep,keep]
}

#######################################################
#######################################################

# Load required modules
# setwd("C:/Users/brand/Desktop/streamlined_ERC_pfam")
args = commandArgs(trailingOnly=TRUE)

# Create a "trees object" containing a single tree for each gene (or domain in your case).
#	Input file format is "name, tab, Newick tree string, newline".
#	Species names must be consistent throughout.

trees=readTrees("ERC/yeast_tree_file_complete.txt", rearrange=T, computePaths=T)

# save(trees, file="yeast_tree_file.RData")	# Save as R data for rapid loading next time.


# Compute full ERC matrix
#	Testing may be done with the "maxDo" argument to limit the number of values to compute.
#	species.threshold sets a minimum number of species to be SHARED between a pair of genes, otherwise no ERC value is computed.
#	species.list allows user to study a subset of species
#load("yeast_tree_file.RData")
#erc=correlateTreesAll(trees,  usePaths=T, maxDo=10)   # Stops after 10 ERC values are calculated for troubleshooting

erc=correlateTreesAll(trees,  usePaths=T, species.threshold=3)

# erc=correlateTreesAll(trees,  usePaths=T, species.threshold=3, species.list= args)
order(erc)

# erc=removeEmpties(erc)

save(erc, file="ERC/yeast_erc_matrix_fixed--test");		# Save ERC matrix for R

write.table( t(erc), "ERC/erc_yeast_fixed--test.tsv", sep="\t", row.names=F, col.names=T , quote=F)		# writes tab-delimited table for other applications

# small function to clean out empty rows and columns from the ERC matrix resulting when genes have an insufficient number of species.



# Code to study a correlation between a specified pair of genes
# This will give you enough data to compute the ERC value for a specified pair instead of the full genome-wide matrix.
#	The "res" is a projection object containing the original branch values "l1 l2", the normalized branch-specific rates "e1 e2", the normalization vector "nv", etc...
#res=correlateTreesProj("gene1", "gene2", trees, plot=T, usePaths = T)
#	This plots the normalized values of gene1 versus gene2 to see the rates that led to the ERC value.
#plotWtext(res$e1, res$e2, labels=res$names, cex=0.7)
