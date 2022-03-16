# SARS-CoV-2 amplicon sequencing scripts
- Repo contains code to:
	1. Assemble SC2 amplicon sequences using the V2, V3, VarSkip, and Midnight primer kits
		- Mainly useful for running on an HPC using Linux and SLURM workload manager
	2. Perform QC and amplicon drop-out analysis
		- Useful for getting a snap-shot of sequencing quality and possible amplicon drop-out issues

- Repo is currently in progress

## pipeline [in progress]
- No script added yet

## consensus-stats.py
- Script to perform a QC of consensus sequences and check for long stretches of N's or possible drop out issues
- Quick usage: `./run-stats.py -c midnight-consensus-seqs.fa -l 1000 --kit midnight -o ./mid-consensus-qc.png`
