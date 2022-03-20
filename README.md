# SARS-CoV-2 amplicon sequencing scripts
- Repo contains code to:
	1. Assemble SC2 amplicon sequences using the V2, V3, VarSkip, and Midnight primer kits
		- Mainly useful for running on an HPC using Linux and SLURM workload manager
	2. Perform QC and amplicon drop-out analysis
		- Useful for getting a snap-shot of sequencing quality and possible amplicon drop-out issues

- Repo is currently in progress

## ancov-pipeline.sh
- Script to run the artic-ncov2019 pipeline on data generated using commonly used kits on an HPC
	- This will need to be run with the `sayf-artic` environment activated; can be done using `conda env create -f sayf-artic.yml`

### Quick usage
- ``


## consensus-stats.py
- Script to perform a QC of consensus sequences and check for long stretches of N's or possible drop out issues

### Quick usage
- `./run-stats.py -c midnight-consensus-seqs.fa -l 1000 --kit midnight -o ./mid-consensus-qc.png`

## amp-drop-analysis.py
- Script to check for amplicon drop-outs given a primer scheme, primer reference, and FASTA/multi-FASTA

### Quick usage
- `./amp-drop-analysis.py -b ~/Documents/TestSpace/midnight/V1/midnight.scheme.bed -r ~/Documents/TestSpace/midnight/V1/midnight.reference.fasta -c ~/Documents/SC2Genomes/Broad-consensus.mfa`
- Example output from above command: ![image](https://user-images.githubusercontent.com/50063071/159175779-064525f4-e0e0-4d60-9da0-6174ea393e6e.png)
