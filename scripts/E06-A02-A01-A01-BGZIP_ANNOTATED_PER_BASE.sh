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
	THREADS=$8
	SAMPLE_SHEET=$9
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=${10}

# Use bgzip to compress the padded and annotated refseq select cds + omim per base depth report

START_PER_BASE_BGZIP=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec ${ALIGNMENT_CONTAINER} bgzip" \
			CMD=${CMD}" -@ ${THREADS}" \
			CMD=${CMD}" -c" \
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt" \
		CMD=${CMD}" >| ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt.gz" \
		CMD=${CMD}" &&" \
		CMD=${CMD}" singularity exec ${ALIGNMENT_CONTAINER} tabix" \
			CMD=${CMD}" -c \#" \
			CMD=${CMD}" -s 1" \
			CMD=${CMD}" -b 2" \
			CMD=${CMD}" -e 2" \
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt.gz"
		CMD=${CMD}" &&" \
			CMD=${CMD}" rm -rvf" \
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt"

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

END_PER_BASE_BGZIP=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},H01,REFSEQ_PER_BASE_BGZIP,${HOSTNAME},${START_PER_BASE_BGZIP},${END_PER_BASE_BGZIP} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
