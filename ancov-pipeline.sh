#!/bin/bash

# <Summary> Script to run the whole artic-ncov2019 pipeline

# Call the right conda source
#source /home/sayf/miniconda3/etc/profile.d/conda.sh

# Get the start-time
st=$(date +%s)

# Set up the option parser
bname=$(basename $0) # Pull the filename

# Define the help function
function help {
	echo " "
	echo "------------------------------"
	echo "| ${bname} 'help' |"
	echo "------------------------------"

	echo "NOTE: Please activate 'sayf-artic' environment prior to running script"
	echo "${bname} -n <read-name> -f <path/to/fastq_pass> -r <path/to/reference.fa> -p <int specifying the kit>"
	#-s <path/to/sequencing-summary>
	#echo "-o: directory name to store the analysis in (default: ./analysis)" # Stored as DIR
	echo "-n: name of the sequencing run"
	echo "-f: Directory containg FASTQ files for each barcode"
	#echo "-b: CSV containing barcodes,samplename"
	#echo "-s: path to sequencing summary (usually in the same directory that fastq/fast5_pass are also in)"
	#echo "-r: Multi-fasta containing reference sequences for MSA and tree-building"
	echo "-p: Primer kits to use. (0: v3, 1:v4, 2: VarSkip, 3: Midnight)"
	echo " "
}

# Create the option parser
clear # Clear the output
	# Put in the default arguments
ODIR="./"

optstring=":hf:n:r:p:"
while getopts ${optstring} opt
do
	case ${opt} in
		h)
			help
			exit 0;;
		#o)
			#ODIR="${OPTARG}"
			##echo "Directory location to deposit the analysis into: ${ODIR}";;
		n)
			RNAME="${OPTARG}"
			echo "Run name '${RNAME}'";;
		f)
			FDIR="${OPTARG}"
			echo "Directory containg FASTQ files for each barcode: ${FDIR}";;
		r)
			REF="${OPTARG}"
			echo "Reference file: ${REF}";;
		#s)
			#SEQSUMM="${OPTARG}"
			#echo "Sequencing summary: ${SEQSUMM}";;
		p)
			SWITCH="${OPTARG}"
			case ${SWITCH} in
				0)
					echo "Artic v3 primer kit chosen"
					;;
				1)
					echo "Artic v4 primer kit chosen"
					;;
				2)
					echo "VarSkip primer kit chosen"
					;;
				3)
					echo "Midnight primer kit chosen"
					;;
			esac
			;;

		?)
			echo ""
			echo "Invalid input: -${OPTARG}"
			help
			exit 1;;
		:)
			echo "Must provide input"
			exit 1;;
	esac
done
echo ""

# Ensure that required arguments are present (RNAME, REF, BCODES)
EXT=0
if [[ -z $RNAME ]]
then
	echo "Please specify a 'Run name' using the '-n' option. Please use ./${bname} '-h' for options."
	EXT=1
fi

#if [[ -z $REF ]]
#then
	#echo "Please specify a reference file using the '-r' option. Please use ./${bname} '-h' for options."
	#EXT=1
#fi

if [[ -z $SWITCH ]]
then
	echo "Please specify which primer kits were used with the '-p' option. Please use ./${bname} '-h' for options."
	EXT=1
fi

if [[ $EXT -eq 1 ]]
then
	exit 1
fi

# Make the directory if it does not already exist
if [[ ! -d ${ODIR} ]]
then
	mkdir ${ODIR}
fi

# Store the current directory
PWD=$(pwd)

# Create the tmp file to figure out when all the jobs are done
TMP=".job_ids.tmp"
echo "TMP Jobs file" > ${PWD}/${TMP}

#TOT_JOBS=$(ls ${PWD}/${FDIR} | grep -c barcode)
TOT_JOBS=0

# Iterate through the files and submit the jobs in parallel
NTHREADS=8
for b in ${PWD}/${FDIR}/barcode*
do
	# Check for at least 20x coverage (~600 BP reads and ~ 30 KB genoms = ~1000 reads for 20x coverage); else pass
	if [[ $(zcat $b/*fastq.gz | head -n 1001 | wc -l) -gt 1000 ]] && [[ $b != "unclassified" ]] # Make sure there is more than one just the barcode; sample attached to barcode
	then
		TOT_JOBS=$((TOT_JOBS + 1))
		bcode=$(python -c "b='$b'; print(b.split('/')[-1])")
		#sname=$(grep $bcode ${PWD}/${BCODES} | awk -F "," '{print $2}')
		sbatch -c 4 --wrap "/home/sayf/pipelines/artic-ncov/single-pipe.sh ${PWD}/${FDIR} ${bcode} ${PWD}/${ODIR} ${RNAME} ${NTHREADS} ${PWD}/${SEQSUMM} ${PWD}/${TMP} ${SWITCH}"
	fi
done

# Cat all of the consensus sequences together
	# - First, we have to wait until all the commands finish running
SLEEP_TIME=0
while [[ $(grep -c barcode ${PWD}/${TMP}) != ${TOT_JOBS} ]]
do
	WAIT_MINS=$(python -c "print(round(${SLEEP_TIME}/60, 2))")
	printf "Waiting ... ${WAIT_MINS} mins passed\r"
	SLEEP_TIME=$((SLEEP_TIME + 30))
	sleep 30
done

echo ""
echo "Concatting all consesnsus sequences"
cat ${PWD}/*.consensus.fasta > ${PWD}/${RNAME}-consensus-seqs.fa &&
cat /home/sayf/pipelines/artic-ncov/references/*f*a ${PWD}/${RNAME}-consensus-seqs.fa > ${PWD}/${RNAME}-consensus-var-msa-seqs.fa &&
echo "Done!"

# Run pangolin
echo ""
echo "Running pangolin"
#/home/sayf/miniconda3/condabin/conda activate pangolin
#srun -c 2 pangolin ${PWD}/${RNAME}-consensus-seqs.fa -t 4 \
	#--outfile ${PWD}/${RNAME}-lineage-report.csv
	#--usher \
#conda deactivate
#echo "Done!"

cat slurm*out > "ancov19-$(date +%m-%d-%y).slurm.log"
rm slurm*.out -f

# Make everything neater; store the extra files in 'IntermediateFiles' directory
mkdir IntermediateFiles
#mv *fastq?* IntermediateFiles
mv *preconsensus* IntermediateFiles
#mv *muscle.in* IntermediateFiles
mv *.vcf* IntermediateFiles
mv *trimmed*bam* IntermediateFiles
mv *primersitereport* IntermediateFiles
#mv *primertrimmed* IntermediateFiles
#mv *sorted.bam? IntermediateFiles
mv *coverage_mask* IntermediateFiles
#mv *muscle.out* IntermediateFiles
mv *minion.log* IntermediateFiles
mv *barcode*fastq IntermediateFiles


if [[ ! -d AlignReports ]]
then
	mkdir AlignReports
fi
mv *alignreport* AlignReports

echo ""
echo "Concatting all consesnsus sequences"
cat ${PWD}/*.consensus.fasta > ${PWD}/${RNAME}-consensus-seqs.fa &&
cat /home/sayf/pipelines/artic-ncov/references/*f*a ${PWD}/${RNAME}-consensus-seqs.fa > ${PWD}/${RNAME}-consensus-var-msa-seqs.fa &&
echo "All consensus sequences have been concattenated!"
#rm *consensus.fasta

# Run MUSCLE prior to IQTree2
echo "Running FAMSA for multiple sequence alignment"
srun -c 12 /home/sayf/utils/famsa -t 24 ${PWD}/${RNAME}-consensus-var-msa-seqs.fa ${PWD}/${RNAME}-msa-seqs.afa &&
#mafft --retree 1 --maxiterate 10 --globalpair --thread 32 ${PWD}/${RNAME}-consensus-var-msa-seqs.fa > ${PWD}/${RNAME}-msa-seqs.fa &&
#/usr/bin/muscle -in ${PWD}/${RNAME}-consensus-seqs.fa -out ${PWD}/${RNAME}-msa-seqs.fa -maxiters 4 -diags

if [[ $? != 0 ]]
then
	echo "MSA failed"
	exit 1
fi
echo "MSA done!"

# Run IQTree2 so we can conduct a phylogenetic comparison of the isolates we have
	# Dow we want to also add in other known lineages to see exactly where we fall, or can pangolin do that?
echo ""
echo "Running IQTree2 for tree inference"
srun -c 2 /home/sayf/utils/iqtree -s ${PWD}/${RNAME}-msa-seqs.afa -m MFP -nt 4 -quiet \
	-o "Wuhan-Hu-1_MN908947.3" -redo
echo "Tree inference is done!"


rm "./core" -f

# Check the coverage stats
case ${SWITCH} in
	0) # V3
		SCHEME=/home/sayf/utils/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed
		kit="v3"
		;;
	1) # V4
		SCHEME=/home/sayf/utils/artic-ncov2019/primer_schemes/nCoV-2019/V4/SARS-CoV-2.scheme.bed
		kit="v4"
		;;
	2) # VarSkip
		SCHEME=/home/sayf/utils/artic-ncov2019/primer_schemes/NEB_VarSkip/V1a/NEB_VarSkip.scheme.bed
		kit="varskip"
		;;
	3) # midnight
		SCHEME=/home/sayf/utils/artic-ncov2019/primer_schemes/midnight/V1/midnight.scheme.bed
		kit="midnight"
		;;
	esac


if [[ ! -d ConsensusSeqs ]]
then
	mkdir ConsensusSeqs
fi
mv *consensus.fasta ConsensusSeqs
#rm -f ${PWD}/${RNAME}-consensus-var-msa-seqs.fa

if [[ ! -d Images ]]
then
	mkdir Images
fi
mv *png Images


if [[ ! -d BamFiles ]]
then
	mkdir BamFiles
fi
mv *bam* BamFiles

rm *hdf -f

echo "Calculating statistics"
srun /home/sayf/SayfGit/SC2/run-stats.py \
	-d ./BamFiles/ -s ./ConsensusSeqs/ -n ${RNAME} \
	-b ${SCHEME} --msa ${PWD}/${RNAME}-msa-seqs.afa &&

srun /home/sayf/SayfGit/SC2/consensus-stats.py -c ${PWD}/${RNAME}-consensus-seqs.fa -l 800 --kit ${kit} -o "${RNAME}-consensus-stats.png" &&
echo "Coverage statistics have been calculated!"

et=$(date +%s)
elapsed=$(python -c "elapsed = round(($et - $st)/60,2); print(elapsed)")
echo "The pipeline completed for ${TOT_JOBS} samples in ${elapsed} mins"
