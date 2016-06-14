#A helper function to parse raw MAF files
#This function reads the provided gene name, sample id, and mutation type and then generates
#a table of counts for the number of mutations each sample has in each gene

#We note that additional parsing of MAF files is necessary in a few cases:
# 1. some different sample IDs refer to the same patient
# 2. some gene names are poorly annotated
#This simple parser does not address these issues

parseMaf = function(protNames, mafFile, mutTypes=c("Missense_Mutation", "Nonsense_Mutation")) {
	ids = c()
	tCount = c()
	for(line in readLines(mafFile)){
	    #find the protein name
	    split = strsplit(line, "\t")[[1]]
	    protInd = match(split[1], protNames) 

	     #check if protein name is in list
	    if(is.na(protInd)){next}
	    
	    #check mutation type
	    if(!any(sapply(mutTypes, function(x) grepl(x,line)))){next} 

	    #find sample id and add to list
	    id = split[grepl("TCGA", split)][1] 
	    idInd = match(id, ids)
	    if(is.na(idInd)){
	        ids = c(ids, id)
	        idInd = length(ids)
	        tCount = rbind(tCount, rep(0, length(protNames)))
	    }

	    #add to count for protein/sample pair
	    tCount[idInd, protInd] = tCount[idInd, protInd]+1
	}
	tCount
}
