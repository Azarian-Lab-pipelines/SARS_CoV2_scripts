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
	1. Prints out alignments for partial primer matches if they exist or notifies you if no match exists for a primer
	2. Generate a plot of the percentage of samples an amplicon does not match
		- Example plot: ![image](https://user-images.githubusercontent.com/50063071/159177273-2465fb87-5367-4d8f-8a65-4772b2658e48.png =100x100) 

### Quick usage
- `./amp-drop-analysis.py -b ~/Documents/TestSpace/midnight/V1/midnight.scheme.bed -r ~/Documents/TestSpace/midnight/V1/midnight.reference.fasta -c ~/Documents/SC2Genomes/Broad-consensus.mfa`
	- Example output from above command: ![image](https://user-images.githubusercontent.com/50063071/159175779-064525f4-e0e0-4d60-9da0-6174ea393e6e.png)

- `./amp-drop-analysis.py -b ~/Documents/TestSpace/midnight/V1/midnight.scheme.bed -r ~/Documents/TestSpace/midnight/V1/midnight.reference.fasta -c ~/Documents/SC2Genomes/Broad-consensus.mfa --text-print | less`
	- The `--text-print` option removes colors from the print outs and allows for neatly saving output in a text file
	- Example output from above command: ![image](https://user-images.githubusercontent.com/50063071/159175853-a1fba359-8418-4562-87cd-c8f55b7f6421.png)
