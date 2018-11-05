## load the list of R data frames produced by evec2SimilarityDf.py
load("L1.gzip")
## connect to Bioconductor and install the necessary packages
source('https://bioconductor.org/biocLite.R')
biocLite('WGCNA')
library(WGCNA)
biocLite('flashClust')
library(flashClust)
## allow multithreading
allowWGCNAThreads(nThreads=12)
## set the value of beta for the similarity matrix to 5
powers <- rep(5,length(L1))
## initialize the list of outputs
l <- vector("list", length(L1))
for(j in c(1:length(L1))){
    sim=L1[[j]]
    GeneNames = colnames(sim)
    adj = adjacency.fromSimilarity(abs(as.matrix(sim)),power=powers[j])
    TOM=TOMsimilarity(adjMat = adj)
    dissTOM = 1-TOM
    #while(ii<11 && sum(dynMod==1)>1){
    geneTree = flashClust(as.dist(dissTOM),method='average')
    dynMod = cutreeDynamic(dendro=geneTree,minClusterSize = 5,method='tree',deepSplit = TRUE)
    l[[j]] = list(geneTree,dynMod,GeneNames[dynMod==1])
    myL <-sapply(GeneNames[dynMod==1],function(nm) sub('[.][.]','',strsplit(nm,"....",fixed=TRUE)[[1]][2]))
    write(myL,sprintf("eigengene_%i_stage.txt",j),sep="\n")
}


