#! /usr/bin/env python

#########################################################
# File: tabulateSumaclust.py                            #
# Author: Shyam Gopalakrishnan                          #
# Date: 15th February 2016                              #
# Description: This script takes the output fasta from  #
# sumaclust and converts it to a tab separated text file#
# with the counts for each sample and each otu.         #
# Optionally, it adds another column with the best blast#
# results for the OTU sequence.                         #
#########################################################

import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate an text file table from sumaclust OTU fasta.")
    parser.add_argument("-i", "--inFasta", dest="fasta", type=str, metavar="InputFasta", required=True, help='Input file (from sumaclust output)')
    parser.add_argument("-blast", "--blast", action="store_true", help="Get best matches from NCBI using blast")
    parser.add_argument("-bold", "--bold", action="store_true", help="Get best matches using BOLD database")
    parser.add_argument("-s", "--scale", dest="scale", type=int, metavar="NumReads", required=False, help='Number of reads to scale each sample to', default=0)
    parser.add_argument("-o", "--out", dest="out", type=str, metavar="outfile", required=False, help="Output text filename", default="")
    args = parser.parse_args()
    if args.out == "":
        args.out = "SampleVsOTUs.txt"
    if args.blast and args.bold:
        print ("You can use only one of -blast and -bold.")
        print ("Choosing the blast option.")
        args.bold = False
    if args.bold:
        print ("The -bold option is not implemented yet.")
        args.bold = False


infile = open(args.fasta)
otunum = {}
countMatrix = {}
otuSeq = {}
seq = ""
samples = []
clusterCenter = False
for line in infile:
    line = line.strip()
    if line[0] == '>':
        toks = line[1:].split()
        if toks[5] == "cluster_center=True;": ## Cluster center, so save its sequence
            clusterCenter = True
        else:
            clusterCenter = False
        ## if seq was not empty then we were at cluster center,
        ## so save the seq
        if (seq != ""):
            otuSeq[otunum[clustname]] = seq
            #print otunum[clustname], clustname, seq
            seq = ""
        (sample, count, clustname) = (toks[0].split(":")[0], int(toks[1].split('=')[1][0:-1]), toks[3].split('=')[1][0:-1])
        if clustname not in otunum:
            otunum[clustname] = "OTU"+str(len(otunum)+1)
        if otunum[clustname] not in countMatrix: ## First occurrrence 
            countMatrix[otunum[clustname]] = {}
        if sample not in samples:
            samples.append(sample)
        if sample in countMatrix[otunum[clustname]]:
            countMatrix[otunum[clustname]][sample] += count
        else:
            countMatrix[otunum[clustname]][sample] = count
    else:
        ## only sequence, so add it and continue. :)
        ## this accounts for multiline sequences in the
        ## fasta
        if clusterCenter:
            seq += line
## Final otu sequence added to list
if seq != "":
    otuSeq[otunum[clustname]] = seq
            
infile.close()
print ("Read data for", len(samples), "samples and", len(otunum), "OTUs.")

if args.scale > 0:
    for sample in samples:
        sampleTotal = 0
        for otu in range(1,len(otunum)+1):
            if sample in countMatrix["OTU"+str(otu)]:
                sampleTotal += countMatrix["OTU"+str(otu)][sample]
        for otu in range(1,len(otunum)+1):
            if sample in countMatrix["OTU"+str(otu)]:
                countMatrix["OTU"+str(otu)][sample] = int(countMatrix["OTU"+str(otu)][sample]*args.scale*1.0/sampleTotal)
    print ("Scaled all samples to", args.scale, "reads.")

outfile = open(args.out, "w")
if args.blast:
    blastfile = open(args.out+".blast.txt", "w")
outfile.write("OTU\t")
outfile.write("\t".join(samples))
outfile.write("\tSeq\n")
for index in range(1,len(otunum)+1):
    curOtu = "OTU"+str(index)
    curCounts = countMatrix[curOtu]
    outline = [curOtu]
    for sample in samples:
        if sample in curCounts:
            outline.append(str(curCounts[sample]))
        else:
            outline.append("0");
    outline.append(otuSeq[curOtu])
    outfile.write("\t".join(outline))
    outfile.write("\n")
    if args.blast:
        blastfile.write(">"+curOtu+"\n")
        blastfile.write(otuSeq[curOtu]+"\n")
outfile.close()
if args.blast:
    blastfile.close()
    print ("Blast can be run using the blast input fasta file:", args.out+".blast.txt")
    
