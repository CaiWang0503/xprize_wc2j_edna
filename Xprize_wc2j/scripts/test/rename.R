#!/usr/bin/Rscript

args<-commandArgs(T)
  fileswarm <- args[1] #Must be an output file from swarm
  min_reads <- args[2]
  outputfolder<-args[3]
  library(Biostrings)
  library(tidyverse)
  fileswarm="/Users/caiwang/Desktop/test/bigtree_6libraris/test.fasta"
  # Read cluster list database
  message("Reading swarm database...")
  db <- readDNAStringSet(fileswarm)
  db <- as.character(db)
  df <- data.frame(SequenceID = names(db), Sequence = db)
  rm(db)
  df=df%>%separate(SequenceID, c("x1","x2","x3","xx4","5","x6","x7","x8","x9"),sep=" ")
  df$x1=paste0("OTU",c(1:length(df$x1)))
  df=df%>%separate(x9, c("barcode","sample"),sep="=")
  df=df%>%separate(x3, c("size","sizenum"),sep="=")
  df$OTUname = paste0(df$sample,":",df$x1,";","size=",df$sizenum)
  message("writing output")
  seq<-apply(df, 1, function(x){paste0(">",x[13],"\n",x[12])})
  write.table(seq,file=paste0(outputfolder,"filtered_seqkit_rename.fasta"),row.names = F,col.names = F,quote=F)
  
 