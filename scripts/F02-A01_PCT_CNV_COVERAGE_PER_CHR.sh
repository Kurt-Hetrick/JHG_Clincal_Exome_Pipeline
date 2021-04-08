# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on

	set

	echo

# INPUT VARIABLES

	CNV_CONTAINER=$1
	ALIGNMENT_CONTAINER=$2
	CORE_PATH=$3

	PROJECT=$4
	FAMILY=$5
	SM_TAG=$6

	CNV_CALL_PCT_BED=$7 # example command to create bed file
		# sed '1d' RefGene.Coding.Clean.Format.19Feb2018_gene_exomeDepth_noY.bed |sort -k1,1 -k2,2n -k3,3n -u|bedtools merge -i - >RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region.bed
		# BEDFILE=/mnt/research/statgen/DDL/CNV_pipeline_DDL/data/RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region2.bed
		# I think I know what this is...I also think it should be created from the original bed file.

# Create header in output

	echo -e "CHROMOSOME\tPERCENT_COVERED_BY_CNVS" \
	>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".call.pct.by.chr.txt"

# function to create an array for each chromosome of its length in the exomeDepth input bed file that has been merged

	CREATE_CHR_LENGTH_ARRAY ()
	{
		CHR_LENGTH_ARRAY=(`awk 'BEGIN {OFS="\t"} $1=="'$CHROMOSOME'" {print $1,($3-$2)}' \
			$CNV_CALL_PCT_BED \
			| sort -k 2,2n \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				-g 1 sum 2`)

			CHROMOSOME_BED=${CHR_LENGTH_ARRAY[0]}

			CHROMOSOME_LENGTH=${CHR_LENGTH_ARRAY[1]}
	}

# function to merge the exomeDepth output per chromosome
# intersect with the exomeDepth input that has been merged
# calculate the total length of the the two above and divide by th total length possible

	CALCULATE_PCT_LENGTH ()
	{
		awk -v CHROMOSOME_LENGTH="$CHROMOSOME_LENGTH" 'BEGIN {OFS="\t"}  $1=="'$CHROMOSOME_BED'" {print $1,$2,$3}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".exomeDepth.bed" \
		| sort -k 2,2n \
		| singularity exec $CNV_CONTAINER \
			bedtools \
			merge -i - \
		| singularity exec $CNV_CONTAINER \
			bedtools \
			intersect \
			-a - \
			-b $CNV_CALL_PCT_BED \
		| awk '{SUM += $3-$2} END {printf "chr" "'$CHROMOSOME_BED'" "\t" "%.1f\n", 100*SUM/"'$CHROMOSOME_LENGTH'"}'
	}

# loop through each unique chromosome in the exomeDepth bed file that has been merged.
# keeping only those chromosomes that start with a number or X or Y
# run both function above
# do a natural sort of karyotype chromosomes

	for CHROMOSOME in $(cut -f 1 $CNV_CALL_PCT_BED | sort | uniq | egrep "^[0-9]|^X|^Y")
	do
		CREATE_CHR_LENGTH_ARRAY
		CALCULATE_PCT_LENGTH
	done \
	| sort -k 1,1V \
	>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".call.pct.by.chr.txt"
