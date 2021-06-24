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

	GATK_3_7_0_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	REF_GENOME=$6
	CODING_BED=$7
		CODING_BED_NAME=$(basename ${CODING_BED} .bed)
		CODING_MD5=$(md5sum ${CODING_BED} | cut -c 1-7)
	PADDING_LENGTH=$8
	GENE_LIST=$9
	SAMPLE_SHEET=${10}
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=${11}

### --Depth of Coverage ALL REFSEQ CODING EXONS PLUS "X" BP FLANKS--
### the padding comes from the GUI

START_DOC_CODING=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec ${GATK_3_7_0_CONTAINER} java -jar"
			CMD=${CMD}" /usr/GenomeAnalysisTK.jar"
		CMD=${CMD}" -T DepthOfCoverage"
			CMD=${CMD}" -R ${REF_GENOME}"
			CMD=${CMD}" -I ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}.bam"
			CMD=${CMD}" -geneList:REFSEQ ${GENE_LIST}"
			CMD=${CMD}" -L ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.bed"
			CMD=${CMD}" -mmq 20"
			CMD=${CMD}" -mbq 10"
			CMD=${CMD}" --outputFormat csv"
			CMD=${CMD}" -ct 10"
			CMD=${CMD}" -ct 15"
			CMD=${CMD}" -ct 20"
			CMD=${CMD}" -ct 30"
			CMD=${CMD}" -ct 50"
		CMD=${CMD}" -o ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD"
		##### MOVE AND RENAME OUTPUTS
		# FIX DELIMETER IN THE HEADER AND ONLY KEEP THE FIRST 2 COLUMNS for the every base report.
		# this report is fed into the annotate per base reports (filtered and unfiltered)
		# DELETE SOURCE FILE
		CMD=${CMD}" &&"
			CMD=${CMD}" sed 's/\t/,/g' ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD"
			CMD=${CMD}" | cut -d \",\" -f 1,2"
			CMD=${CMD}" >| ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.EveryBase.csv" 
			CMD=${CMD}" &&"
				CMD=${CMD}" rm -rvf ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD"
		# COUNT OF BASES AT "X" LEVEL OF COVERAGE
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_cumulative_coverage_counts"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_cumulative_coverage_counts.csv"
		# FRACTION OF BASES AT "X" LEVEL OF COVERAGE
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_cumulative_coverage_proportions"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_cumulative_coverage_proportions.csv"
		# SUMMARY COVERAGE STATISTICS FOR EACH GENE WHERE PADDED CODING BED FILE OVERLAPS GENE LIST
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_gene_summary"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_gene_summary.csv"
		# COUNT OF INTERVALS COVERAGE BY AT LEAST "X" LEVEL OF COVERAGE
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_interval_statistics"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_interval_statistics.csv"
		# SUMMARY STATISTICS FOR EACH INTERVAL IN PADDED CODING BED FILE.
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_interval_summary"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_interval_summary.csv"
		# NOT SURE WHAT THIS IS AT THE MOMENT
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_statistics" \
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_statistics.csv"
		# SAMPLE COVERAGE SUMMARY STATISTICS FOR SAMPLE IN PADDED CODING BED FILE.
		CMD=${CMD}" &&"
			CMD=${CMD}" mv -v ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_summary"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/$FAMILY/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.sample_summary.csv"

	# write command line to file and execute the command line

		echo ${CMD} >> ${CORE_PATH}/${PROJECT}/COMMAND_LINES/${SM_TAG}_command_lines.txt
		echo >> ${CORE_PATH}/${PROJECT}/COMMAND_LINES/${SM_TAG}_command_lines.txt
		echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			if [ "${SCRIPT_STATUS}" -ne 0 ]
				then
					echo ${SM_TAG} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
					>> ${CORE_PATH}/${PROJECT}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
					exit ${SCRIPT_STATUS}
			fi

END_DOC_CODING=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},E01,DOC_CODING,${HOSTNAME},${START_DOC_CODING},${END_DOC_CODING} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
