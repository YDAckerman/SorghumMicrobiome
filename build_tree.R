library(seqinr)
library(ape)

adata <- read.alignment("~/Documents/PurdomGroup/data/otu.fasta",
                        format = "fasta")

## taken from http://a-little-book-of-r-for-bioinformatics.readthedocs.io/
##                 en/latest/src/chapter5.html#building-an-unrooted-
##                 phylogenetic-tree-for-protein-sequences

mymat  <- as.matrix.alignment(adata)
alignment <- ape::as.alignment(mymat)
alignmentbin <- as.DNAbin(alignment)
mydist <- dist.dna(alignmentbin)
mytree <- nj(mydist)
mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.


## the Susan Holmes way:
## source("https://bioconductor.org/biocLite.R")
## biocLite("DECIPHER")
library(DECIPHER)
library(phangorn)
alignment <- AlignSeqs(DNAStringSet(unlist(adata$seq)), anchor = NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order

## plot the tree
plot(treeNJ)

save(alignment, treeNJ, file = "phyloTreeData.rda")
