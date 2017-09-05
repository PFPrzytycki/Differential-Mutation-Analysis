#A helper function to parse raw MAF files
#This function reads the provided gene name, sample id, and mutation type and then generates
#a table of counts for the number of mutations each sample has in each gene

#We note that additional parsing of MAF files is necessary in a few cases:
# 1. some different sample IDs refer to the same patient
# 2. some gene names are poorly annotated
#This simple parser does not address these issues

parseMaf = function(protNames, mafFile, usePackages, mutTypes=c("Missense_Mutation", "Nonsense_Mutation")) {
    
  tryCatch({

    if(!usePackages){stop("usePackages")}

    #Read maf file skipping comment lines
    mut = fread(mafFile, skip="Hugo_Symbol", header=TRUE, fill = TRUE)
    
    #Trim proteins and mutation types
    mut = mut[mut$Hugo_Symbol %in% protNames & mut$Variant_Classification %in% mutTypes,]
    
    #If ids are TCGA barcodes, trim them for uniqueness
    if(grepl("TCGA", mut$Tumor_Sample_Barcode[1])){
      mut$Tumor_Sample_Barcode = sapply(mut$Tumor_Sample_Barcode, function(x) paste(strsplit(x, "-")[[1]][1:4], collapse="-"))
    }
    ids = unique(mut$Tumor_Sample_Barcode)

    #Count mutations per patients per gene
    tCount = count(mut, vars=c("Hugo_Symbol","Tumor_Sample_Barcode"))
    
    #Convert counts to matrix
    tCount =  spMatrix(length(ids), length(protNames),
                       match(tCount$Tumor_Sample_Barcode, ids), match(tCount$Hugo_Symbol, protNames),
                       tCount$freq)

    tCount

  }, error = function(err) {
    print("Switching to slow read, either flag usePackages is set to FALSE or could not find one of the following columns: Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode")
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
    return(tCount)
  })

}
