#!/bin/bash

# <Summary> Script to run the whole artic-ncov2019 pipeline for each barcode in a parallel manner

##### Store the variables #####
FDIR=${1} # Full path for fastq_pass
BCODE=${2} # Barcode name
ODIR=${3} # Full path to output_DIR
RNAME=${4} # Run name to save the files with
NTHREADS=${5} # Number of threads
#SNAME=${6} # Sample name
F5DIR=$(echo ${FDIR} | sed 's/q/5/g')
SEQSUMM=${6}
TMP=${7}
SNAME=$(echo $BCODE | sed 's/barcode/sample/g')

# Option for which primer kit to use
SWITCH=${8} # 0: v3 -- 1: v4 -- 2: VarSkip -- 3: Midnight

# Get the start-time
st=$(date +%s)

##### Run the pipeline #####
module purge
module load guppy-cpu
#################### UNSURE ####################
#source /home/sayf/miniconda3/etc/profile.d/conda.sh
#/home/sayf/miniconda3/condabin/conda activate artic-ncov2019

# Define the different paramters as different functions

# 1) De-multiplex the data, one barcode at a time?
	# - Run strict de-multiplexing to avoid the chimeras that arise from ONT

#echo ""
#echo "Demultiplexing ${BCODE}"
#guppy_barcoder -t ${NTHREADS} --require_barcodes_both_ends -i ${FDIR}/${BCODE} -s ${ODIR} --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" && \
	 #grep -c "^@" ${ODIR}${BCODE}/*fastq | awk -v code=${f} -F ":" '{total += $2}END{print "Total number of reads in ", code":", total}' || \
	 #echo "No reads passed for ${BCODE}" &&

# The output will be multiple files; all the fastq's are de-multiplexed and they are then batched into multiple files
	# - Each batch fastq file will contain 4000


# Use the appropriate pipeline based on the primers
case ${SWITCH} in
	0) 
		echo "V3 primer kit used"
		# Perform the size-based read-filtering
		artic guppyplex --skip-quality-check --min-length 400 --max-length 700 \
			--directory ${FDIR}/${BCODE} --prefix ${RNAME} &&

		# Run the pipeline using v4 primers
		artic minion --normalise 200 --threads ${NTHREADS} \
			--read-file ${RNAME}_${BCODE}.fastq --scheme-directory /home/sayf/utils/artic-ncov2019/primer_schemes \
			--fast5-directory ${F5DIR}/${BCODE} --sequencing-summary ${SEQSUMM} \
			nCoV-2019/V3 ${SNAME} &&
			echo "Done!"

		;;

	1) 
		echo "V4 primer kit used"
		# Perform the size-based read-filtering
		artic guppyplex --skip-quality-check --min-length 400 --max-length 700 \
			--directory ${FDIR}/${BCODE} --prefix ${RNAME} &&

		# Run the pipeline using v4 primers
		artic minion --normalise 200 --threads ${NTHREADS} \
			--read-file ${RNAME}_${BCODE}.fastq --scheme-directory /home/sayf/utils/artic-ncov2019/primer_schemes \
			--fast5-directory ${F5DIR}/${BCODE} --sequencing-summary ${SEQSUMM} \
			nCoV-2019/V4 ${SNAME} &&
		echo "Done!"

		;;

	2) 
		echo "VarSkip primer kit used"
		# Perform the size-based read-filtering
			# We skip the quality check because we are only using the fastq_pass reads already
		artic guppyplex --skip-quality-check --min-length 500 --max-length 900 \
			--directory ${FDIR}/${BCODE} --prefix ${RNAME} &&

		# Perform the Minion pipeline using the VarSkip primers
		artic minion --normalise 200 --threads ${NTHREADS} \
			--read-file ${RNAME}_${BCODE}.fastq --scheme-directory /home/sayf/utils/artic-ncov2019/primer_schemes \
			--scheme-version 1a NEB_VarSkip ${SNAME} --skip-nanopolish --medaka --medaka-model r941_min_fast_g303 &&
		echo "Done!"

		;;
	3)
		echo "Midnight primer kit used"
		# Perform the size-based read-filtering
			# We skip the quality check because we are only using the fastq_pass reads already
		artic guppyplex --skip-quality-check --min-length 1000 --max-length 1500 \
			--directory ${FDIR}/${BCODE} --prefix ${RNAME} &&

		# Perform the Minion pipeline using the VarSkip primers
		artic minion --normalise 200 --threads ${NTHREADS} \
			--read-file ${RNAME}_${BCODE}.fastq --scheme-directory /home/sayf/utils/artic-ncov2019/primer_schemes \
			--scheme-version 1 midnight ${SNAME} --skip-nanopolish --medaka --medaka-model r941_min_fast_g303 &&
			echo "Done!"
		;;

	?)
		echo "Invalid input: ${SWITCH}"
		;;

	"")
		echo "Must provide input"
		;;
esac


# Perform the Minion pipeline using the VarSkip primers

# Get the end-time and print the elapsed time
et=$(date +%s)
elapsed=$(python -c "elapsed = round(($et - $st)/60,2); print(elapsed)")

#conda deactivate
echo "Pipeline completed in ${elapsed} minutes for ${BCODE}"
echo "${BCODE} done" >> ${TMP}
