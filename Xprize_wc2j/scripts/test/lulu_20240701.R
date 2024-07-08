## ------------------------------------------------------------------------------
# set some general options
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) # if you want to pass arguments to R from the shell (bash) command line
print(args)

otutablefilename <- args[1]

## ------------------------------------------------------------------------------
suppressMessages(library(tidyverse))
# devtools::install_github("tobiasgf/lulu", force = TRUE)
library(lulu)

## ---- eval=FALSE, include=FALSE------------------------------------------------
## ------------------------------------------------------------------------------
otutable_raw <- read_tsv(otutablefilename) %>% column_to_rownames(var = "OTU")
otutable=otutable_raw
otutable <- otutable %>% dplyr::select(-Seq)  

matchlist <- read.table("matchlist.txt", header=FALSE,
                        as.is=TRUE, stringsAsFactors=FALSE)

## ------------------------------------------------------------------------------
curated_result <- lulu(otutable, matchlist, minimum_match = 84) # default minimum_match = 84% similarity, but such a low value only works if there are many samples, so that otu co-occurrence gives a strong signal that two otus are parent and child.

curated_table <- rownames_to_column(curated_result$curated_table, 
                                    var = "OTU"
                                    )

write_tsv(curated_table, file = "table_97_lulu84.txt")

###get lulu sequence######
otutable_raw$OTU=rownames(otutable_raw)
newotutab=left_join(curated_table,otutable_raw,by="OTU") %>% select(OTU,Seq)
seq<-apply(newotutab, 1, function(x){paste0(">",x[1],"\n",x[2])})
write.table(seq,file="OTUs_sumaclust_lulu84_forblast.fasta",row.names = F,quote=F)
#need to delete x from xxx.fasta by hand.
write.csv(newotutab, file = "./data/10outof1/otutable_97_10outof1_LULU.csv",row.names = TRUE, col.names = F,quote = FALSE)