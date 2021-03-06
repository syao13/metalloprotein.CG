#!/mlab/data/software/R-3.2.1-F22/bin/Rscript

##!/usr/bin/Rscript

############################################################################################
##
##   Written by Robert Flight, 07/20/2016
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
##
###########################################################################################

## ------------------------------------------------------------------------
library(atom2seq)
library(Biobase)
library(Biostrings)

## ----testAlign-----------------------------------------------------------
# need to add a test for non-standards existing as insertion character.
atom <- AAStringSet(c(rmFrontAddBack="MDDTEKMSMKLT", addFrontRmBack="TSMDDTEKMSMK", insertion="MDDAATEKMSMKL", deletion="SMTEKMSMKL", indel="SMTEKMSMAAKL", nonAAInsert="SMDD--TEKMSMKL", nonAAInsertRmFront="MDD--TEKMSMKL"))

seq <- AAStringSet(c(rep("SMDDTEKMSMKL", 5), nonAAInsert="SMDDTEKMSMKL", nonAAInsertRmFront="SMDDTEKMSMKL"))

refNum <- list(rmFrontAddBack=c(2,3,4,5,6,7,8,9,10,11,12,0),
               addFrontRmBack=c(0,1,2,3,4,5,6,7,8,9,10,11),
               insertion=c(2,3,4,0,0,5,6,7,8,9,10,11,12),
               deletion=c(1,2,5,6,7,8,9,10,11,12),
               indel=c(1,2,5,6,7,8,9,10,0,0,11,12),
               nonAAInsert=c(1,2,3,4,0,0,5,6,7,8,9,10,11,12),
               nonAAInsertRmFront=c(2,3,4,0,0,5,6,7,8,9,10,11,12))
names(seq) <- names(atom)
useNames <- names(atom) # do a small sample for testing

## ----genAlignTest--------------------------------------------------------
outTest <- lapply(useNames, function(x){genAlign(atom[[x]], seq[[x]])})
outNum <- lapply(outTest, function(x){x$atomNum})
names(outNum) <- useNames
all.equal(outNum, refNum)

## ----loadData------------------------------------------------------------
#useDir <- "/mlab/data/rmflight/Documents/projects/work/sen/coordination_families"
#setwd("../output_allMetal/")
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

## ----read_sequences----------------------------------------
library(Biostrings)
seq <- readAAStringSet("seqs.SEQRES.shell.txt", use.names=TRUE)
atom <- readAAStringSet("seqs.ATOM.shell.txt", use.names=TRUE)

## ----sanityCheck-------------------------------------------
splitDot <- function(inNames){
  allSplit <- strsplit(inNames, ".", fixed=TRUE)
  outNames <- sapply(allSplit, function(x){x[[1]]})
}
seqNames <- splitDot(names(seq))
atomNames <- splitDot(names(atom))

all.equal(seqNames, atomNames)

seqSort <- sort(names(seq))
atomSort <- sort(names(atom))

all.equal(seqSort, atomSort)

## ----runAlign------------------------------------------------
useNames <- sample(names(atom), 100)
outAlign <- mclapply(useNames, function(x){genAlign(atom[[x]], seq[[x]])}, mc.cores=6)

## ----runAllAlign---------------------------------------------
useNames <- names(atom)
outNum <- mclapply(useNames, function(x){
  tmpRes <- genAlign(atom[[x]], seq[[x]])
  tmpRes$atomNum
}, mc.cores=6) ## This takes about 20 min
names(outNum) <- useNames

save(outNum, file="atom2seq.RData")



