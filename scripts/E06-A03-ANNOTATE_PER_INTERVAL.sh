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
# redirecting stderr/stdout to file as a log.

	set

	echo

# INPUT VARIABLES

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	CODING_BED=$6
		CODING_BED_NAME=$(basename ${CODING_BED} .bed)
		CODING_MD5=$(md5sum ${CODING_BED} | cut -c 1-7)
	PADDING_LENGTH=$7

START_PER_INTERVAL=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# convert the per coding interval report in DEPTH_OF_COVERAGE/CODING_PADDED into a bed file with the depth
# intersect with the annotated coding padded bed file

	awk 'NR>1' ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_interval_summary.csv \
		| sed 's/:/\t/g; s/,/\t/g' \
		| awk 'BEGIN {OFS="\t"} \
			$2!~"-" \
			{print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14} \
			$2~"-" \
			{print $0}' \
		| awk 'BEGIN {OFS="\t"} \
			{gsub (/-/,"\t"); \
			print $1,$2,$3,$4,$5,$8,$9,$10,$13,$14,$15}' \
		| singularity exec $ALIGNMENT_CONTAINER bedtools \
			intersect \
			-a ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD-ANNOTATED.bed \
			-b - \
			-wb \
		| awk '{OFS=","} \
			BEGIN {print "#CHROM" "," "START" "," "END" "," "GENE" "," "TRANSCRIPT" "," "EXON" "," "STRAND" \
			"," "TOTAL_CVG" "," "AVG_CVG" "," "Q1_CVG" "," "MEDIAN_CVG" "," "Q3_CVG" "," "PCT_20x" "," \
			"PCT_30x" "," "PCT_50x"} \
			{print $1,$2,$3,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18}' \
	>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.CODING_EXON_SUMMARY.csv

END_PER_INTERVAL=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},F01,REFSEQ_PER_INTERVAL,${HOSTNAME},${START_PER_INTERVAL},${END_PER_INTERVAL} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv
