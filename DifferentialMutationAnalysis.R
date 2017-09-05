###############################################################################
# Differential Mutation Analysis
#
# The only required input is a single MAF file, a sample MAF file for BRCA is
# provided in "Data/BRCA_sample.maf" and can be tested by calling
#
# DifferentialMutationAnalysis("Data/BRCA_sample.maf")
#
# Output is a single two or three column file with protein names, their 
# uEMD scores, and optionally, supporting q-values, named after the input
# file with DiffMut appended to it (e.g. "BRCA_sample-DiffMut.txt")
#
# The code can optionally be run to search for oncogenes or tumor suprressor
# genes separately by passing "onco" or "TSG" as options for geneType
#
# DifferentialMutationAnalysis("Data/BRCA_sample.maf", geneType="onco")
#
# Finally, the code computes supporting q-values for genes. To 
# compute q-values simply pass a value p which determines the numer of
# background distributions to generate (default is 5). Note that this comes at 
# a cost to runtime. The code can be run with no permutations to quickly output
# a list of genes ranked by uEMD score.
# 
# DifferentialMutationAnalysis("Data/BRCA_sample.maf", p=0)
#
# If you have the following three R packages: "data.table", "plyr", and
# "Matrix" you can set the flag usePackages to TRUE to significantly decrease
# file read time
###############################################################################

#helper function to parse MAF file
source("parseMaf.R")

#rank normalize mutaton or variation counts
fastRank = function(x) { 
  x[x!=0] = rank(x[x!=0], ties.method="min")+length(which(x==0)) 
  x/length(x) 
}

#bin counts to build histogram
bins = function(v, p=100) {
  l = length(v)
  nBins = rep(0,p)
  for(val in v){
    nBins[ceiling(val*(p-1))+1] = nBins[ceiling(val*(p-1))+1]+1
  }
  nBins = nBins/l
  nBins
}

#compute unidirectional EMD between mutationas and variations
uEMD = function(tBins, nBins, p=100) {
  sum = 0
  move = 0
  for(i in p:2){
    move = move+tBins[i]-nBins[i]
    sum = sum+max(move,0)
  }
  sum
}

#generate random uEMDs to compute FDRs
generateRandomEMDs = function(tRank, nRankBinned) {
  permRank = t(apply(tRank, 1, sample))
  permRankBinned = apply(permRank, 2, bins)
  randEMDs = sapply(1:dim(nRankBinned)[2], function(x) uEMD(permRankBinned[,x], nRankBinned[,x]))
  randEMDs
}

#compute FDRs based on random uEMDs
computeFDR = function(uEMDscores, randEMDs) {
  FDRs = sapply(uEMDscores, function(x) length(which(randEMDs>=x))/(length(which(uEMDscores>=x))))
  FDRs
}

#compute q-values from FDRs
computeQ = function(FDRs, uEMDscores) {
  qVal = sapply(1:length(FDRs), function(x) min(FDRs[uEMDscores<=uEMDscores[x]]))
  qVal
}

#Main Function for Differential Mutation Analysis
DifferentialMutationAnalysis = function(mafFile, geneType="all", p=5,
                                        outDir = "Output/", 
                                        protFile = "Data/protNames.txt", 
                                        natBinFile = "Data/natDistBinned.txt",
                                        usePackages = FALSE) {

  if(usePackages){
    library("data.table")
    library("plyr")
    library("methods")
    library("Matrix")
  }

  #A list of protein names
  protNames = read.table(protFile, stringsAsFactors=FALSE)$V1

  #load ranked binned natural variation count data
  if(!usePackages){
    nRankBinned = read.table(natBinFile)
  }
  else{
    nRankBinned = fread(natBinFile, data.table=FALSE)
  }

  #determine if we want to find all cancer genes or just oncogenes or TSGs
  if(geneType=="onco"){
    tCount = parseMaf(protNames, mafFile, usePackages, "Missense_Mutation")
  }
  if(geneType=="TSG"){
    tCount = parseMaf(protNames, mafFile, usePackages, "Nonsense_Mutation")
  }
  else{
    tCount = parseMaf(protNames, mafFile, usePackages)
  }
  
  #rank normalize mutations
  tRank = t(apply(tCount, 1, fastRank))

  #bin the rank distribution
  tRankBinned = apply(tRank, 2, bins)

  #compute uEMD scores
  uEMDscore = sapply(1:length(protNames), 
    function(x) uEMD(tRankBinned[,x], nRankBinned[,x]))

  #create output directory if it doesn't exist
  if(!dir.exists(outDir)){ dir.create(outDir) }

  #output only uEMD scores if no q-values are needed (faster run time)
  if(p==0){
    write.table(cbind(protNames, uEMDscore), 
      paste0(outDir, strsplit(basename(mafFile),".maf")[[1]],"-DiffMut.txt"), quote=FALSE, row.names=FALSE)
  }
  else{
    #compute q-values, p determines number of times random uEMDs are generated
    FDRs = rowSums(sapply(1:p, function(x) 
      computeFDR(uEMDscore, generateRandomEMDs(tRank, nRankBinned))))/p
    qVal = computeQ(FDRs, uEMDscore)
    write.table(cbind(protNames, uEMDscore, qVal), 
      paste0(outDir, strsplit(basename(mafFile),".maf")[[1]],"-DiffMut.txt"), quote=FALSE, row.names=FALSE)
  }

}
