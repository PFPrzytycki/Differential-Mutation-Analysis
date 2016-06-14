###############################################################################
# Differential Mutation Analysis
#
# The only required input is a single MAF file, a sample MAF file for BRCA is
# provided in "Data/BRCA_sample.maf" and can be tested by calling
#
# DifferentialMutationAnalysis("Data/BRCA_sample.maf")
#
# Output is a single two column file with protein names and their uEMD scores
# named "uEMDscores.txt"
#
# The code can optionally be run to search for oncogenes or tumor suprressor
# genes separately by passing "onco" or "TSG" as options for geneType
#
# DifferentialMutationAnalysis("Data/BRCA_sample.maf", geneType="onco")
#
# Finally, the code can optionally compute supporting p-values for genes. To 
# compute p-values for the p highest scoring genes, simply pass a value p
#
# DifferentialMutationAnalysis("Data/BRCA_sample.maf", p=100)
#
# This will generate a separate p-value file ordered by uEMD score called
# "pvalues.txt"
#
###############################################################################

#helper function to parse MAF file
source("parseMaf.R")

#rank normalize mutaton or variation counts
fastRank = function(x) { 
  x[x!=0] = rank(x[x!=0], ties.method="min")+length(which(x==0)) 
  x/length(x) }

#bin counts to build histogram
bins = function(v, p=100){
  l = length(v)
  nBins = rep(0,p)
  for(val in v){
    nBins[ceiling(val*(p-1))+1] = nBins[ceiling(val*(p-1))+1]+1
  }
  nBins = nBins/l
  nBins
}

#compute unidirectional EMD between mutationas and variations
uEMD = function(t_v, n_v, p=100){
    tBins = bins(t_v, p)
    nBins = bins(n_v, p)
    sum = 0
    move = 0
    for(i in p:1){
        move = move+tBins[i]-nBins[i]
        sum = sum+max(move,0)
    }
    sum
}

#compute p-value for gene i
pvalsP = function(i, uEMDscore, nRank, p, perms=100) {
  length(which(sapply(1:perms, function(x) 
    uEMD(sample(0:99, p, TRUE, bins(nRank[,i])+1/p^2)/99, nRank[,i]))>=uEMDscore[i]))/perms
}

#Main Function for Differential Mutation Analysis
DifferentialMutationAnalysis = function(mafFile, geneType="all", p=0){
  #A list of protein names
  protNames = read.table("Data/protNames.txt", stringsAsFactors=FALSE)$V1

  #determine if we want to find all cancer genes or just oncogenes or TSGs
  if(geneType=="onco"){
    tCount = parseMaf(protNames, mafFile, "Missense_Mutation")
  }
  if(geneType=="TSG"){
    tCount = parseMaf(protNames, mafFile, "Nonsense_Mutation")
  }
  else{
    tCount = parseMaf(protNames, mafFile)
  }
  #load natural variation count data
  natDistMis = read.table("Data/protMisDists_p3.txt")
  
  #rank normalize mutations and variations
  nRank = t(apply(natDistMis, 2, fastRank))
  tRank = t(apply(tCount, 1, fastRank))

  #compute uEMD scores
  uEMDscore = sapply(1:length(protNames), function(x) uEMD(tRank[,x], nRank[,x]))

  write.table(cbind(protNames, uEMDscore), "uEMDscores.txt", quote=FALSE, row.names=FALSE)

  #optionally compute p-values
  if(p>0){
    pval = sapply(order(-uEMDscore)[1:p], function(i) pvalsP(i, uEMDscore, nRank, dim(tRank)[1]))
    write.table(cbind(protNames=protNames[order(-uEMDscore)[1:p]], pval), "pvals.txt", quote=FALSE, row.names=FALSE)
  }
}
