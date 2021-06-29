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

	GATK_3_7_0_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	REF_GENOME=$5
	SAMPLE_SHEET=$6
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=$7

# APPLY VQSR INDEL MODEL TO THE VCF

START_APPLY_RECALIBRATION_INDEL=`date '+%s'`

	# construct command line

		CMD="singularity exec ${GATK_3_7_0_CONTAINER} java -jar"
			CMD=${CMD}" /usr/GenomeAnalysisTK.jar"
		CMD=${CMD}" --analysis_type ApplyRecalibration"
			CMD=${CMD}" --reference_sequence ${REF_GENOME}"
			CMD=${CMD}" --disable_auto_index_creation_and_locking_when_reading_rods"
			CMD=${CMD}" --mode INDEL"
			CMD=${CMD}" --ts_filter_level 99.9"
			CMD=${CMD}" --input:VCF ${CORE_PATH}/${PROJECT}/TEMP/CONTROLS_PLUS_${FAMILY}.VQSR.SNP.vcf"
			CMD=${CMD}" --recal_file ${CORE_PATH}/${PROJECT}/${FAMILY}/VCF/VQSR/CONTROLS_PLUS_${FAMILY}.HC.INDEL.recal" \
			CMD=${CMD}" --tranches_file ${CORE_PATH}/${PROJECT}/${FAMILY}/VCF/VQSR/CONTROLS_PLUS_${FAMILY}.HC.INDEL.tranches"
		CMD=${CMD}" --out ${CORE_PATH}/${PROJECT}/TEMP/CONTROLS_PLUS_${FAMILY}.VQSR.SNP_INDEL.vcf"

	# write command line to file and execute the command line

		echo ${CMD} >> ${CORE_PATH}/${PROJECT}/COMMAND_LINES/${FAMILY}_command_lines.txt
		echo >> ${CORE_PATH}/${PROJECT}/COMMAND_LINES/${FAMILY}_command_lines.txt
		echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			if [ "${SCRIPT_STATUS}" -ne 0 ]
				then
					echo ${FAMILY} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
					>> ${CORE_PATH}/${PROJECT}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
					exit ${SCRIPT_STATUS}
			fi

END_APPLY_RECALIBRATION_INDEL=`date '+%s'`

# write out timing metrics to file

	echo ${FAMILY}_${PROJECT},J01,APPLY_RECALIBRATION_INDEL,${HOSTNAME},${START_APPLY_RECALIBRATION_INDEL},${END_APPLY_RECALIBRATION_INDEL} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
