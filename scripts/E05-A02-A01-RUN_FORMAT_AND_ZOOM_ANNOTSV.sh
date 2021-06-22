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
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	
	FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT=$6
	ZOOM_LIST=$7
	ZOOM_NAME=$8

	SAMPLE_SHEET=$9
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=${10}

## run R script for exomeDepth

START_FORMAT_AND_ZOOM_ANNOTSV=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec ${CNV_CONTAINER} Rscript" \
			CMD=${CMD}" ${FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT}" \
			CMD=${CMD}" --annotSVtemp ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}.annotSV.temp.tsv" \
			CMD=${CMD}" --zoomfile ${ZOOM_LIST}" \
			CMD=${CMD}" --zoomname ${ZOOM_NAME}" \
			CMD=${CMD}" --smTag ${SM_TAG}" \
			CMD=${CMD}" --outdir ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CNV_OUTPUT"

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

END_FORMAT_AND_ZOOM_ANNOTSV=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},G01,FORMAT_AND_ZOOM_ANNOTSV,${HOSTNAME},${START_FORMAT_AND_ZOOM_ANNOTSV},${END_FORMAT_AND_ZOOM_ANNOTSV} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from samtools bam to cram

	exit ${SCRIPT_STATUS}
