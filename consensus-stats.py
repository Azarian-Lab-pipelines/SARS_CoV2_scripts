#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import sys


#################### Set up the argument parser ####################
import argparse as ap
parser = ap.ArgumentParser(description="""
Script to summarize consensus sequence quality and potential drop-out issues
""", formatter_class=ap.RawTextHelpFormatter)

parser.add_argument("-c", "--consensus", dest="consensus", required=True,
        help="Specify the consensus sequence file")

parser.add_argument("-l", "--length", dest="drop_len", required=True, type=int,
        help="""Specify the minimum length of contiguous N's to consider as potential drop-out.
        We recommend that the region of interest be slightly smaller 
        than the amplicon size to ensure that possible drop-outs are 
        not missed due to overlapping amplicon schemes""")

parser.add_argument("--kit", dest="kit",
        choices=["v3", "v4", "varskip", "midnight"],
        help="Choose which of the provided primer kits you are using")
parser.add_argument("--bed", dest="bed", 
        help="Specify the BED file containing the primer scheme if not choosing one of the built-in kits")
parser.add_argument("--kit-name", dest="kit_name",
        help="Specify the kit name if a BED file is provided")

parser.add_argument("-o", "--output", dest="out", required=True,
        help="Specify the output name to save the image (*.png|tiff|pdf|svg|jpg)")

args = parser.parse_args()

# Make sure necessary arguments are present
if not args.kit and not args.bed:
    parser.error("No primer kit has been chosen and no BED file has been provided. One of the two must be provided. See --help for options")

if args.bed and not args.kit_name:
    parser.error("No kit name provided but is required when a BED file is provided")


#################### Import the gene and amplicon location files ####################
vadr_url = "https://raw.githubusercontent.com/sayfaldeen/SC2/main/wgs-vadr-genes.bed?token=GHSAT0AAAAAABSFAPAII6ROMI3UUHPISW72YR3RICQ"
genes = pd.read_csv(vadr_url, sep="\t",
                   names=["Ref", "Start", "End", "Gene"])

amp_url = "https://raw.githubusercontent.com/sayfaldeen/SC2/main/amplicon-positions.tsv?token=GHSAT0AAAAAABSFAPAIHDVU5DYVTORHCQFUYR3RK6A"
df = pd.read_csv(amp_url, sep="\t")

if args.kit:
    amp_df = df[df.Scheme == args.kit]
    kit = args.kit

else:
    kit = args.kit_name
    bed = pd.read_csv(args.bed, sep="\t",
                 header=None,
                 names=["Ref", "PS", "PE", "Primer", "Pool"])

    bed["Amplicon"] = bed.Primer.apply(lambda x:int(re.search("_\d+_", x).group().strip("_")))

    bed = bed.sort_values(by="Amplicon", ascending=True).reset_index(drop=True)

    res = {}
    n=0
    for a in bed.Amplicon.unique():
        tmp = bed[bed.Amplicon == a]
        res[n] = [a, kit, tmp.iloc[0]["PS"], tmp.iloc[1]["PE"]]
        n += 1

    amp_df = pd.DataFrame.from_dict(res, orient="index",
                                    columns = ["Amplicon", "Scheme", "Start", "End"])


#################### Set up the necessary functions ####################
def FindSeqs(f):
    # Pull the sequences
    seqs = {}
    seq = None
    
    for line in open(f, "r"):
        if line.startswith(">"):
            
            if seq:
                seqs[header] = "".join(seq)
            
            header = line.strip()
            seq = []
            
        else:
            seq.append(line.strip())

    seqs[header] = "".join(seq) # The final sequence will not have a header after it
    
    return seqs


def FindStretches(F, drop_len = 100):
    """
    Function to find stretches of N's
    """
    
    seqs = FindSeqs(F)
    global df_cols 
    df_cols = ["Genes", "Locations", "Lens", f"{kit.capitalize()}Amplicon", "N-percent"]
    
    rdf = pd.DataFrame(columns=df_cols)
    for h, seq in seqs.items():

        # 1) Pull the areas with stretches of 'N'
        matches = dict([(stretch.span(), stretch.group()) for stretch in re.finditer("N{" + str(drop_len) +",}", seq)])

        # 2) Find genes and amplicons with stretches in them
        # Corrected mega-list comps
        gene_stretches = [set([genes[(genes.Start < s) & (genes.End > s)]["Gene"].values[0], genes[(genes.Start < e) & (genes.End > e)]["Gene"].values[0]]) \
            for s,e in [loc for loc in matches.keys()]\
            if len(genes[(genes.Start < s) & (genes.End > s)]["Gene"].values) > 0 and \
            len(genes[(genes.Start < e) & (genes.End > e)]["Gene"].values) > 0 ]

        c_amp = [set([amp_df[(amp_df.Start < s) & (amp_df.End > s)]["Amplicon"].values[0], amp_df[(amp_df.Start < e) & (amp_df.End > e)]["Amplicon"].values[0]]) \
            for s,e in [loc for loc in matches.keys()]\
            if len(amp_df[(amp_df.Start < s) & (amp_df.End > s)]["Amplicon"].values) > 0 and \
            len(amp_df[(amp_df.Start < e) & (amp_df.End > e)]["Amplicon"].values) > 0 ]
        

        # 3) Add the N-content
        n_cont = round(seq.count("N")/len(seq) *100,2)
        
        # 4) Add the location and length into the dataframe
        ndf = pd.DataFrame.from_dict({h: [gene_stretches, 
                    [x for x in matches.keys()],
                  [e-s for s,e in matches.keys()],
                                     c_amp, n_cont]},
                                     orient='index',
                                    columns=df_cols)
        rdf = rdf.append(ndf)

    rdf["Drops"] = rdf.Lens.apply(lambda x:len(x) if len(x) > 0 \
            else 0)
    rdf["Sample"] = rdf.apply(lambda x:x.name.split("/")[0].strip(">"), axis=1)
    rdf.set_index("Sample", inplace=True)

    global drops
    drops = rdf[rdf.Drops > 0].copy()
    drops["CumeLens"] = drops.Lens.apply(lambda x:sum(x))

    print("")
    print(f"{len(drops)}/{len(rdf)} sequences have stretches of greater than {drop_len} contiguous N's")
    
    return rdf

#################### Run the functions to create the DF ####################
rdf = FindStretches(F=args.consensus, drop_len=args.drop_len)

if len(drops) == 0:
    print("Since no samples have any drop-outs no plots will be generated")
    sys.exit(0)


#################### Make and save the plots ####################
# Prep the data vectors for the plots
len_drops = []
drops.Lens.apply(lambda x:len_drops.extend(x))
gene_drops = []
drops.Genes.apply(lambda x:gene_drops.extend(x))

# This will need to be amended to allow choice of amplicon
amp_drops = []
drops[df_cols[3]].apply(lambda x:[amp_drops.extend(list(v)) for v in x])

# Make the plots
plt.figure(figsize=(14,16), 
          facecolor="white")
plt.rcParams.update({'axes.facecolor':'white'})

plt.style.use("ggplot")
plt.suptitle("Descriptive plots of consensus sequence 'N' sites\n\n",
            size=14, weight='bold')

plt.subplot(2,2,1)
sns.histplot(x=len_drops, bins=np.arange(0, 17000, step=500),
            color="lightblue", edgecolor="k")
plt.xticks(np.arange(0, 17000, step=1_000),
          rotation=45, ha='right')
plt.tick_params(axis='x', zorder=2, width=1, size=6)
plt.title(f"Distributon of the lengths of stretches with > {args.drop_len} contiguous N's\n")
plt.xlabel("\nCumulative length of stretches")
plt.ylabel("Count")

plt.subplot(2,2,2)
sns.countplot(x=[" >--< ".join(x) for x in gene_drops],
             color="lightblue", edgecolor="k")
plt.title(f"Count of genes with stretches of > {args.drop_len} contiguous N's\n")
plt.xlabel("\nGene")
plt.xticks(rotation=45, ha='right')
plt.ylabel("Count")

plt.subplot(2,2,3)
sns.histplot(x=drops.Drops.values,
            color="lightblue", edgecolor="k")
plt.title(f"Number of stretches with > {args.drop_len} contiguous N's\n")
plt.xlabel("\nNumber of stretches")
plt.ylabel("Count")

# Add in plot for amplicon drop-out sites
plt.subplot(2,2,4)
sns.countplot(x=amp_drops,
             color="lightblue", edgecolor="k")
plt.title(f"{kit} amplicons with > {args.drop_len} contiguous N's\n")
plt.xlabel(f"{kit} Amplicon")
plt.ylabel("Count")


plt.tight_layout(w_pad=2, h_pad=6)
plt.savefig(args.out, dpi=200)

# Print out that the script has completed
print(f"Image saved as '{args.out}'")
