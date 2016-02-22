#!/usr/bin/Rscript

##############  This program is to get non-redundant set of zinc ID list. ############## 
## Non-redundant set is defined by the shell domain as well as binding ligand combination. 
## If the shell ligands spans over multi-chains, consider them all together.

####### find multi metal within redundant set


options(stringsAsFactors=FALSE)
#setwd("../output_allMetal")
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

library(Biostrings)
aaCodes <- toupper(AMINO_ACID_CODE)

## SEQRES sequences, names of atom2seq are shell headers
bindFasta <- readAAStringSet("seqs.SEQRES.bind.txt", format="fasta")
shellFasta <- readAAStringSet("seqs.SEQRES.shell.txt", format="fasta")
load("atom2seq.RData")

length(bindFasta)
length(shellFasta)
bindHeaders <- names(bindFasta)
shellHeaders <- names(shellFasta)

## Order the headers so that sequences are comparable
getId <- function(header) {
  items<-strsplit(header, '[|]')[[1]]
  paste(items[1], substring(tail(items, n=1), 1,1), sep="." )
}

bindId <- sapply(bindHeaders, getId)
shellId <- sapply(shellHeaders, getId)

shellFasta1 <- shellFasta[shellId %in% bindId]
length(shellFasta1)

shellHeaders1 <- names(shellFasta1)
shellId1 <- sapply(shellHeaders1, getId)


bindFasta <- bindFasta[order(bindId)]
shellFasta <- shellFasta1[order(shellId1)]
bindHeaders <- names(bindFasta)
shellHeaders <- names(shellFasta)


## put together sequences and corresponding header for binding ligands and shell ligand
#sum(bindFasta==shellFasta)
sequences <- sapply(1:length(bindFasta), function(x) as.character(bindFasta[[x]]))
seqMat <- cbind(bindHeaders, shellHeaders, sequences)


########################## obtain shell sequence domains ############################
## given header, return certain info about multichain zinc:
## min num, max num, min lig, max lig, and ordered by chain
headerBoundaryLig <- function (header) {
  items <- strsplit(header, '[|]')[[1]]
  
  ## chain of zinc id, and chain of first atom num of the sequence
  chains <- c(strsplit(items[1], '[.]')[[1]][2], substring(items[length(items)], 1,1)) 
  
  itemSpt <- strsplit(items[2:(length(items)-1)], '[.]')
  resIDs <- unlist(sapply(itemSpt, function(x) { if (x[2] %in% aaCodes) {x[1]} }))
  resChains <- substring(resIDs, 1,1)
  resNums <- substring(resIDs, 2)
  resNames <- unlist(sapply(itemSpt, function(x) { if (x[2] %in% aaCodes) {x[2]} }))
  
  chainRg <-sapply(sort(unique(resChains)), function(x) { 
    range(as.numeric(substring(resIDs[which(resChains == x)], 2))) })
  
  rgLig <- sapply(sort(unique(resChains)), function(x) {
    c(resNames[which(resIDs == paste(x, chainRg[1,x], sep=""))[1]],
      resNames[which(resIDs == paste(x, chainRg[2,x], sep=""))[1]] )} )
  
  rbind(chainRg, rgLig)
}

#headerBoundaryLig("1B55.A.1|A143.HIS.ND1.N|A154.CYS.SG.S|A155.CYS.SG.S|A165.CYS.SG.S|A2 ")
#headerBoundaryLig("1D1T.A.403|B271.HIS.NE2.N|A507.CYS.OXT.O|B2377.NAD.N7A.N|B735.HOH.O.O|A1")   

### get shell sequnce from seqMat rows and mapping info, outNum 
getOneChainDomain <- function (seqMatItem, outNum) {
  header <- seqMatItem[2]
  seq <- seqMatItem[3]
  
  rangeOfLig <- headerBoundaryLig(header)
  
  seqChain <- substring(tail(strsplit(header, '[|]')[[1]], n=1), 1,1)
  chainStart <- as.numeric(substring(tail(strsplit(header, '[|]')[[1]], n=1), 2))
  
  if (! seqChain %in% colnames(rangeOfLig)) { return("error: seq chain not match ligand chain") }
  
  mapping <- outNum[seqMatItem[2]]
  if (length(mapping[[1]]) == 0) {
    return("error: no mapping")
    next
  }
  
  min <- as.numeric(rangeOfLig[1,seqChain])
  max <- as.numeric(rangeOfLig[2,seqChain])
  
  start <- mapping[[1]][min - chainStart + 1]
  end <- mapping[[1]][max - chainStart + 1]
  if (length(start) != 1 || length(end) != 1) {
    return("error, out of mapping range")
    next
  }
  
  seqMinLig <- as.character(substring(seq, start, start))
  seqMaxLig <- as.character(substring(seq, end, end))     
  mapLigMin <- aaCodes[seqMinLig]
  mapLigMax <- aaCodes[seqMaxLig]
  if (is.na(mapLigMin) ||  is.na(mapLigMax) || mapLigMin != rangeOfLig[3,seqChain] || mapLigMax != rangeOfLig[4,seqChain] ) {
    return("error, wrong mapping")
    next
  }  
  
  start5 <- start - 5
  end5 <- end + 5
  if (start5 < 1) {start5 <- 1}
  if (end5 > nchar(seq)) {end5 <- nchar(seq)}
  
  substring(seq, start5, end5)
}


#getOneChainDomain(seqMat[14235,], outNum)
#        sequences
#"LIGQIHDALHLANEC"

## get all metal-binding domains
domains <- apply(seqMat, 1, function(x) getOneChainDomain(x, outNum))
#length(domains)
#length(unique(domains))

err <- domains[which(substring(domains, 1,5) == "error")]
table(err)

#hist(nchar(domains[which(substring(domains, 1,5) != "error")]), 100,
#     xlab="zinc-binding shell domain length", main="The length of zinc binding shell domain")
#table(nchar(domains))
#mean(nchar(domains))
#sd(nchar(domains))

#################### get unique shell-domain-binding-ligand combinations #####################
## Acquire ligand combination
ligandsFromHeader <- function(header) {
  items <- strsplit(header, '[|]')[[1]]
  
  ligSplit <- strsplit(items[2:(length(items)-1)], '[.]')
  ligNames <- sapply(ligSplit, function(x){x[2]})
  
  paste(sort(ligNames), collapse = '.')
}

#ligandsFromHeader("1B55.A.1|A143.HIS.ND1.N|A154.CYS.SG.S|A155.CYS.SG.S|A165.CYS.SG.S|A2 ")
#ligandsFromHeader("1D1T.A.403|B271.HIS.NE2.N|A507.CYS.OXT.O|B2377.NAD.N7A.N|B735.HOH.O.O|A1")

ligands <- sapply(seqMat[,1], function(x) ligandsFromHeader(x))
#length(ligands)

seqLigs <- paste(domains, ligands, sep=".")
#length(unique(seqLigs))

## combine multiple chains, if any
metalids <- sapply(seqMat[,1], function(x) strsplit(x, '[|]')[[1]][1])
#length(metalids)
uid <- unique(metalids)
#length(uid)
## the reason that uid is not 17135 is because when generating fasta file for multichains,
## removed zincs with all non-aa ligands

## concatenate seq-lig combos together for multi chains
## id, the metalic id to be combine
## metalids, all zinc ids of bindFasta or shellFasta, same length as seqLig, and serves as index for seqLig
## seqLig, single seq-lig combos to be combined by zinc ids
mtSeqFromMetalid <- function(id, metalids, seqLig) {
  i <- which(metalids == id)
  if (sum(substring(seqLig[i],1,5) == "error") > 0) {"error"}
  else if (length(i) == 1) {seqLig[i]}
  else {paste(seqLig[i], collapse=".")}
}

multiSeqLigs <- sapply(uid, function(x) mtSeqFromMetalid(x, metalids, seqLigs))
length(multiSeqLigs)

uMultiSeqLigs <- unique(multiSeqLigs)[unique(multiSeqLigs) != "error"]
#print("Total non-redundant set is 6501")
length(uMultiSeqLigs)
## 6501
## uMultiSeqLigs and finalMetalList have the same dimention


######################################################################
## finalMetalList is the non-redundant list
## size is the number of each redundant set
## no resolution filter included
##
## get the best one for each unique seq-lig combo 
## if no x-ray, use latest, 
## if x-ray only, use the min resolution one,
## if both x-ray and NMR, if minRes < 2, use the minRes one, otherwise, use latest one
######################################################################
#metal <- c(rep("Metal", 18114), rep("Mg", 27340), rep("Ca", 19398), rep("Fe", 14353), rep("Na", 10419))

alldata <- NULL
for (metal in c("Zn", "Mg", "Ca", "Fe", "Na")) {
  setwd(paste("/mlab/data/sen/projects/metal/go_ec_analysis/", metal, "_results/root_afterRF/", sep=""))
  load("rf.results.RData")
  rawdata <- read.table("r.allLig.txt", header=FALSE)

  id <- rawdata[,1]
  orderid <- order(id)
  data <- rawdata[orderid,]

  alldata <- rbind(alldata, cbind(data, rep(metal, dim(data)[1]), prediction.all))
}


#colnames(alldata) <- c("metalID", "method", "year", "resolution", "angleCombo",
#                       "ligandCombo", "bondlengthCombo", "biStatusCombo", "bfactorCombo",
#                       "biLigs", "chainIDCombo", "residueCombo", "atomCombo", "extra", "metal")

rawdata <- alldata


multiMetalList <- NULL
size <- NULL
for (i in 1:length(uMultiSeqLigs)) {
  combo <- uMultiSeqLigs[i]
  if ( substring(combo, 1,5) == "error") {next}
  
  ind <- which(multiSeqLigs==combo)
  size <- c(size, length(ind))
  metalId <- uid[ind] ## uid is required, which is unique metal ids and is one-to-one matching to multiSeqLigs
  
  if (length(metalId)==1) { 
    next
  }
  
  if (FALSE %in% (metalId %in% rawdata[,1])) {next}  

  mets <- sapply(metalId, function(x) rawdata[which(rawdata[,1] == x), 15])
  pred <- sapply(metalId, function(x) rawdata[which(rawdata[,1] == x), 16])
  if (length(table(mets)) == 1) {
    next
  } else {
    multiMetalList <- rbind(multiMetalList, cbind(metalId, mets, pred, rep(i, length(mets))))
  }
}

colnames(multiMetalList) <- c("metalID", "metal", "comp1norm2", "redundantGroup")
rownames(multiMetalList) <- NULL

write.table(multiMetalList, file="/mlab/data/sen/projects/metal/go_ec_analysis/allMetal_results/beforeNR/multiMetalList.txt")

################ below is some characterization of the final list ################
## size of each non-redundant set
#length(finalMetalList)
#hist(size[size!=1], 100, xlab = "size of the redundant sets", ylab="Count")
#sum(size==1)
#table(size)
#mean(size)
#sd(size)

## mothod, year, resolution of the non-redundant set
#mthYrRes <- t(sapply(finalMetalList, function(x) idMthYrRes[which(idMthYrRes[,1] == x), 2:4]))  

## resolution
#data<- rawdata[rawdata[,1] %in% finalMetalList, ]
#sum(data[,24]<3) ##5573 
#hist(as.numeric(data[,24]),100)
#mean(as.numeric(data[,24]))
#sd(as.numeric(data[,24]))

## minimum angles
#minAngle <- apply(data[,2:6], 1, min)
#hist(as.numeric(minAngle,100))


####### check some of the redundant list 
## uid and multiSeqLigs have one-to-one match
## size, finalMetallist, and uMultiSeqLigs have one-to-one match
## metalids and seqMat rows have one-to-one match
#table(size)
#which(size == 25) 
#i <- 1193 ## choose one number from above results
#size[i]
#finalMetalList[i]
#uMultiSeqLigs[i]
#uid[which(multiSeqLigs == uMultiSeqLigs[i])]

#idMthYrRes[(idMthYrRes[,1] %in% uid[which(multiSeqLigs == uMultiSeqLigs[i])]),]
#seqMat[(metalids %in% uid[which(multiSeqLigs == uMultiSeqLigs[i])]),]


###### prepare for Robert
#length(multiSeqLigs)

data<- rawdata[rawdata[,1] %in% finalMetalList, ]
res3metalList <- data[data[,4]<3, 1]

#print("With resolution greater than 3: 6199")
length(res3metalList)

metalListSeqPos <- NULL
for (i in res3metalList) {
  idx <- which(metalids == i)
  
  for (j in idx) {
    seq <- seqMat[j, 3]
    shellH <- seqMat[j,1]
    
    items <- strsplit(shellH, '[|]')[[1]]
    chain <- substring(tail(items,1), 1, 1)
    chainStart <- as.numeric(substring(tail(items,1), 2))
    
    ligs <- NULL
    poss <- NULL
    for (k in 2:(length(items)-1)) {
      if (substring(items[k],1,1) != chain) {next}
      if (! strsplit(items[k], '[.]')[[1]][2] %in% aaCodes) {next}
          
      mapping <- outNum[seqMat[j,2]]
     
      ligNum <- as.numeric(substring(strsplit(items[k], '[.]')[[1]][1], 2))
      ligPos <- mapping[[1]][ligNum - chainStart + 1]
      
      seqLig1 <- as.character(substring(seq, ligPos, ligPos))
      seqLig3 <- aaCodes[seqLig1]
      if (is.na(seqLig3) || seqLig3 != strsplit(items[k], '[.]')[[1]][2] ) {
        #print(shellH)
        #print(items[k])
        next
      }  
      
    ligs <- c(ligs, seqLig3)
    poss <- c(poss, ligPos)
    }
    
    metalListSeqPos <- rbind(metalListSeqPos, 
                          c(i, length(idx), as.character(which(idx==j)), paste(ligs, collapse="."), 
                            paste(poss, collapse="."), as.character(seq)))
  }
}

#dim(metalListSeqPos)
colnames(metalListSeqPos) <- c("metal id", "total chains", "chain counts", "ligand combo", "ligand position combo", "sequence")

#metalListSeqPos[1:5,]

write.table(metalListSeqPos, file="nonredundant.list.txt", sep="\t")











