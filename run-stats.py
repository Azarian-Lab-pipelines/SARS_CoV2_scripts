#!/usr/bin/env python

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

import pysam

# Ignore warning
import warnings
warnings.filterwarnings("ignore")

#################### Set up parser ####################
import argparse as ap
parser = ap.ArgumentParser()
parser.add_argument("-d", "--directory", dest="d",
        help="Specify directory containing the sorted BAM files")
parser.add_argument("--pattern", dest="pattern", default="sample[0-9]*.sorted.bam$",
        help="Specify the naming pattern of the sorted BAM files (default: 'sample[0-9]*.sorted.bam$')")
parser.add_argument("-n", "--name", dest="name", default="",
        help="Specify the name of the run")
parser.add_argument("-b", "--bed", dest="bed",
        help="Specify a BED file for the primers used")
parser.add_argument("--msa", dest="msa_file",
        help="Specify the file containing the MSA")

args = parser.parse_args()

if not args.d:
    parser.error("Please indicate a directory using the '-d' option. See --help for options.")
if not args.bed:
    parser.error("Please indicate a BED file using the '-b' option. See --help for options.")

#################### Find the files ####################
pattern = args.pattern.replace("*", "+").replace("$", "\Z")

if args.d[-1] != "/":
    args.d = args.d + "/"

files = [args.d + f for f in os.listdir(args.d) if re.match(pattern, f)]

#################### Pull data ####################
# Combine the coverages for all of the files
df = pd.DataFrame()
plt.figure(figsize=(15,12))

coverages = {}
for color, f in zip(sns.color_palette("husl", len(files)), files):

    # Pull the coverages for each file
    df[f.split("/")[-1].split(".sorted.bam")[0]] = [int(line.split("\t")[2]) for line \
                               in pysam.depth("-a", f).strip().split("\n")]

    # Calculate the PCT ID for each read
    #pysam.index(f)
    plt.style.use("bmh")
    sam = pysam.AlignmentFile(f, "rb")
    pct_ids = [round((read.query_alignment_length/read.query_length)*100,2) \
           for read in sam.fetch()]

    sns.kdeplot(pct_ids, color=color, 
                label=f.split("/")[-1].split(".")[0], 
               alpha=0.7) 
    
#################### Create images ####################
print(f"Generating images...")

plt.title("Distribution of percents of query mapped to reference\n")
plt.xlabel("\nPercent of query mapped to reference")
plt.ylabel("Density\n")
plt.legend(bbox_to_anchor=(1.1,1, 0.2, 0),
      ncol=2, frameon=True, edgecolor='k', facecolor='white',
      fontsize="large")

plt.savefig(f"{args.name}-read-pctid.png", dpi=200,
           bbox_inches="tight")
print("PCT ID plot done")

#################### Plot the coverage by position first ####################
plt.style.use("bmh")
plt.figure(figsize=(15,12))

for color, column in zip(sns.color_palette("husl", len(files)), df.columns):
    sns.lineplot(x=df.index, y=df[column],
                label=column, color=color, alpha=0.8)
    
plt.legend(bbox_to_anchor=(1.1,1, 0.2, 0), fontsize='large',
          ncol=2, frameon=True, edgecolor='k', facecolor='white')

plt.title("Coverage for every position on the genome\n", size=12)
plt.ylabel("Coverage\n")
plt.xlabel("\nPosition")

plt.savefig(f"{args.name}-position-coverage.png", dpi=200,
           bbox_inches="tight")
print("Positonal coverage plot done")

# Pull in the bed file
bed = pd.read_csv(args.bed, sep="\t",
                 header=None,
                 names=["Ref", "Start", "End", "Name",
                       "Pool", "Strand"])

# Sort by the amplicon number
bed["Amplicon"] = bed.Name.apply(lambda x:int(x.split("_")[1]))
bed.sort_values(by="Amplicon", inplace=True)

bed.Amplicon = bed.apply(lambda x:str(x["Amplicon"]) + "_" + x["Name"].split("_")[-1],
                         axis=1)


# Find the regions covered by each amplicon
regs = {}

for i in bed.index:
    if "LEFT" in bed.iloc[i]["Name"]:
        regs[bed.iloc[i]["Amplicon"].split("_")[0]] = (int(bed.iloc[i]["Start"]), 
                                 bed.iloc[i+1]["End"])
        
# Create a dataframe containing the average coverage for regions covered by an amplicon for each sample
adf = pd.DataFrame()

for a in regs:
    adf = adf.append(pd.DataFrame(df.iloc[regs[a][0]:regs[a][1]].median(),
            columns=[a]).T)

try:    
    adf["Amplicon"] = adf.apply(lambda x:int(x.name),
                            axis=1)
except:
    print("ERROR!")
adf.set_index("Amplicon", inplace=True)


#################### Plot the average coverage for the amplicons for each sample as a bar ####################
plt.style.use("bmh")
plt.figure(figsize=(15,12))

for color, column in zip(sns.color_palette("husl", len(files)), adf.columns):
    sns.barplot(x=sorted(adf.index), y=np.log10(adf[column]),
                    label=column, alpha=0.8, linewidth=2,
                    edgecolor=color, facecolor="white")
   
plt.legend(bbox_to_anchor=(1.05,1, 0.18, 0),
          ncol=2, frameon=True, edgecolor='k',facecolor='white',
          fontsize="medium")

plt.title("Amplicon coverage\n", size=12)
plt.ylabel(r"Median coverage ($log_{10}$)"+"\n", size=10)
plt.xlabel("\nAmplicon", size=10)


plt.xticks([x for x in sorted(adf.index)], rotation=45, ha='center',
        size=8)
plt.yticks(size=8)

plt.savefig(f"{args.name}-amplicon-coverage-bar.png", dpi=200,
           bbox_inches="tight")


#################### Plot the average coverage for the amplicons for each sample as a line ####################
plt.style.use("bmh")
plt.figure(figsize=(15,12))

for color, column in zip(sns.color_palette("husl", len(files)), adf.columns):
    sns.lineplot(x=sorted(adf.index), y=np.log10(adf[column].replace(0,1e-1)),
                label=column, alpha=0.8, linewidth=2,
                color=color)
    
plt.legend(bbox_to_anchor=(1.05,1, 0.18, 0),
          ncol=2, frameon=True, edgecolor='k',facecolor='white',
          fontsize="medium")

plt.title("Amplicon coverage\n", size=12)
plt.ylabel(r"Median coverage ($log_{10}$)"+"\n", size=10)
plt.xlabel("\nAmplicon", size=10)


plt.xticks([x for x in sorted(adf.index)], rotation=45, ha='center',
        size=8)
plt.yticks(size=8)

plt.savefig(f"{args.name}-amplicon-coverage.png", dpi=200,
           bbox_inches="tight")


#################### Plot the median coverage per sample ####################
med_covs = df.median().sort_values()
plt.figure(figsize=(15,12))
plt.style.use("bmh")

sns.barplot(x=med_covs.index, y=med_covs.values, palette="husl",
           edgecolor='k', linewidth=1, alpha=0.7)
ymin,ymax = plt.ylim()
x=0
for s,m in dict(med_covs).items():
    if m > ymax/10:
        plt.text(x=x, y=m*0.9, s=m, ha='center',
                bbox=dict(facecolor="white", boxstyle='round,pad=0.5'),
                fontsize=8)
    else:
        plt.text(x=x, y=ymax/40, s=f"{m}", ha='center',
                bbox=dict(facecolor="white", boxstyle='round,pad=0.5'), 
                fontsize=8)
    x+=1

plt.title("Median coverage per sample\n")
plt.ylabel("Median coverage\n")

plt.xticks(rotation=45, ha="right")

plt.savefig(f"{args.name}-median-coverage.png", dpi=200,
           bbox_inches="tight")
print("Median coverage plot done!\n")


#################### Pull the percent of N's in the consensus sequences ####################
def FindSeqs(f, correct=False):
    # Pull the sequences
    seqs = {}
    seq = None
    
    for line in open(f, "r"):
        if line.startswith(">"):
            
            if seq:
                seqs[header] = "".join(seq)
            
            if correct:
                header = line.strip().replace("_1", "")
            else:
                header = line.strip()

            seq = []
            
        else:
            seq.append(line.strip())

    seqs[header] = "".join(seq) # The final sequence will not have a header after it
    
    return seqs

cons_patt = "sample[0-9]+.consensus.fasta\Z"
consensus_files = [f for files in os.listdir("./") if re.match(cons_patt, f)]

if args.msa_file:
    # Visualize the consensus sequences
    cdf = pd.DataFrame.from_dict(FindSeqs(args.msa_file), orient="index",
                      columns=["Seq"])

    plt.style.use("ggplot")
    plt.figure(figsize=(15,12))

    y = 0
    for i in [x for x in cdf.index if "variant" not in x]:
        seq = cdf.loc[i, "Seq"]
        slen = len(seq)
        
        npos = np.asarray([pos for pos,nuc in enumerate(seq) if nuc.lower() == "n"])
        gpos = np.asarray([pos for pos,nuc in enumerate(seq) if nuc == "-"])
        tpos = np.asarray([n for n in range(0, len(seq)+1)])

        y += 0.1
        sns.scatterplot(x=npos, y=y,
                    color='red', label="N", zorder=2,
                       marker="|", size=0.1, edgecolor='red',
                        legend=False)
        try:
            sns.scatterplot(x=gpos, y=y,
                    color='cyan', label="Gap", zorder=2,
                           marker="|", size=0.1, edgecolor="cyan",
                            legend=False)
        except:
            pass

        sns.lineplot(x=tpos, y=y,
                    color='k', label="Genome", zorder=1,
                       linestyle="-", legend=False)


        plt.xlabel("\nPosition")
        
        plt.text(x=slen/2, y=y+(0.05),
                s=i.replace(">", ""), ha="center", size=6,
                bbox=dict(boxstyle='round', facecolor='white', alpha=1))
        

    # plt.xlim(0, slen+1)
    plt.ylim(0.1, y+0.5)

    plt.yticks([])
    plt.xticks(np.arange(0, slen+100, step=1_000),
              rotation=45, ha='right', va='top')

    plt.title("Visualization of consensus sequences\n")

    # Fix the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    nlabs = []
    nhands = []

    for h,l in zip(handles, labels):
        if l not in nlabs:
            nlabs.append(l)
            nhands.append(h)
            
    plt.legend(nhands, nlabs,
              facecolor="white", frameon=True, edgecolor="k",
              ncol=3, fontsize="medium")

    plt.grid(color='white', linewidth=1)

    plt.savefig(f"{args.name}-consensus-msa-plot.png", dpi=200,
           bbox_inches="tight")


#################### Plot N content ####################
def FindNCont(f):
    # Pull the sequences
    reads = []
    for line in open(f, "r"):
        if line.startswith(">"):
            header = line.strip().split("/")[0].replace(">", "")
        else:
            reads.append(line.strip())
            
    seq = "".join(reads)
    return (header, round((seq.count("N")/len(seq)) *100,3))


cons_patt = "sample[0-9]+.consensus.fasta\Z"
try:
    n_cont = dict([FindNCont(f) for f in os.listdir("./") if re.match(cons_patt, f)])

    plt.figure(figsize=(15,12))
    plt.style.use("default")
    plt.style.use("bmh")
    plt.bar(x=n_cont.keys(), height=n_cont.values(), color='red',
           edgecolor="k")

    plt.xticks(rotation=45, ha='right')

    plt.title("N content per sample\n")
    plt.ylabel("Percent N content\n")

except ZeroDivisionError:
    plt.figure(figsize=(15,12))

plt.savefig(f"{args.name}-n-content-plot.png", dpi=200,
        bbox_inches="tight")
#################### Save the images into one PDF ####################
image_files = [f"{args.name}-position-coverage.png",
                f"{args.name}-amplicon-coverage-bar.png",
                f"{args.name}-amplicon-coverage.png",
                f"{args.name}-median-coverage.png",
                f"{args.name}-read-pctid.png",
                f"{args.name}-consensus-msa-plot.png",
                f"{args.name}-n-content-plot.png"]
                

import PIL as pil

images = [pil.Image.open(img).convert("RGB") for img in image_files]

images[0].save(f"{args.name}-coverage-stats.pdf", save_all=True, 
               append_images = images[1:])

print()
#print(f"Genome position coverage graph saved as '{args.name}-position-coverage.png'")
#print(f"Amplicon coverage graph saved as '{args.name}-amplicon-coverage.png'")
#print(f"Median coverage per sample saved as '{args.name}-median-coverage.png'")
#print(f"Distribution of read alignment percentage storaed as '{args.name}-read-pctid.png'")
print(f"Coverage statistics saved as '{args.name}-coverage-stats.pdf'")
print()

# Write out the coverage stats as a table for each sample
pd.DataFrame(df.median(axis=0).T, 
            columns=["Median Coverage"]).to_csv(f"{args.name}-table.tsv", sep="\t")
