#!/usr/bin/env python3

#################### Import the modules ####################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import os
import re

import pysam
from io import StringIO

#################### Import the argument parser #################### nvm; not needed 
#import argparse as ap
#parser = ap.ArgumentParser(description="""
#Script to perform statisitcs on contig assembly.
#Generates plots for:
    #1) Contig lengths
    #2) Genome coverage
#Then converts the SAM files of the reads that mapped to the SC2 genome > FASTA for polishing
#""", formatter_class=ap.RawTextHelpFormatter)
#
#parser.add_argument("")

#################### Set up the functions ####################
def SamToFA(sam):
    """
    The name of the function says it all... oh and sam is a SAM file
    Used to turn SAM files (from contig alignment to SC2 reference) into FASTA.
    These contigs will then be used for polishing and final assembly
    """
    
    with open(sam.replace("sam", "fasta"), "w") as fa:
        fa.write(pysam.fasta(sam))

def PlotCoverages(f):
    """
    Function to create coverage plots from the sorted SAM files
    """
    
    depths = pd.read_csv(
        StringIO(pysam.depth(f, "-aa")), 
        sep="\t", header=None, index_col=1)\
    .drop(0, axis=1)

    # Condense the line-plot to bins
    s = 1 # Start
    dd = {} # dict to hold results

    bin_step = 100 # Size of the bins
    for e in np.arange(bin_step+1, depths.index.max()+1, step=bin_step):
        dd[s] = np.log10(np.median(depths.loc[s:e, 2].values)) # Take the log_10 of the median
        s = e # Move to the next bin
        
    # Plot the coverage
    plt.plot(dd.keys(), dd.values(),
            label=re.search("barcode\d+", f).group(),
            alpha=0.6)
    
    plt.title("Genome Coverage")
    plt.ylabel(r'Log$_{10}$ Coverage')
    plt.xlabel("Genome Position")
    plt.yticks(np.arange(1, max(dd.values()), step=1))
    
    plt.legend()

def PlotContigLengths(label, contig_file):
    """
    Function to plot distribution of contig lengths
    """
    
    headers = dict([(l.split()[0], l.strip().split()[-1]) \
                 for l in open(contig_file, "r") if l.startswith(">")])
    contig_lengths = [int(x.replace("len=", "").strip()) for x in headers.values()]
    
    plt.hist(contig_lengths, label=label, bins=50,
            alpha=0.6)
    
    plt.xlim(0,2000)
    
    plt.title("Contig Lengths")
    plt.xlabel("Contig Length")
    plt.ylabel("Count")
    
    plt.legend()

#################### Run the functions ####################
# os.walk = dirpath, dirnames, filenames

plt.figure(figsize=(16,8), facecolor="white")
plt.subplot(1,2,1)
# Find the contig files
contig_files = {}
for dirpath, dirname, fname in os.walk("./"):
    if "Assemblies" in dirpath:
        sname = dirpath.split("/")[-1].replace("Assemblies","")
        fname = f"{dirpath}/final.contigs.fa"
        contig_files[sname] = fname

# Plot the contig length plots
for label, cf in contig_files.items():
    PlotContigLengths(label, cf)
    
plt.subplot(1,2,2)
# Plot the genome coverages
files = [f for f in os.listdir(".") if re.search("\A.*sorted.*sam\Z", f)]
for f in files:
    PlotCoverages(f)
    SamToFA(f)
plt.savefig("SC2-assembly-coverage-plots.png", dpi=300)
