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

	ANNOVAR_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	ANNOVAR_DATABASE_FILE=$6
	ANNOVAR_REF_BUILD=$7
	ANNOVAR_INFO_FIELD_KEYS=$8
		ANNOVAR_INFO_FIELD_KEYS_SPACED=$(echo ${ANNOVAR_INFO_FIELD_KEYS} | sed 's/,/ /g')
	ANNOVAR_HEADER_MAPPINGS=$9
		ANNOVAR_HEADER_MAPPINGS_SPACED=$(echo ${ANNOVAR_HEADER_MAPPINGS} | sed 's/,/ /g')
	ANNOVAR_VCF_COLUMNS=${10}
		ANNOVAR_VCF_COLUMNS_SPACED=$(echo ${ANNOVAR_VCF_COLUMNS} | sed 's/,/ /g')
	THREADS=${11}
	SAMPLE_SHEET=${12}
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=${13}

## ANNOTATE VARIANT ONLY CFTR REGION VCF WITH GENE/TRANSCRIPT WITH ANNOVAR

START_ANNOVAR=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec ${ANNOVAR_CONTAINER} python"
			CMD=${CMD}" /annovar_wrangler/annovar_wrangler.py"
			CMD=${CMD}" --vcf_input_path ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}.VARIANT_SITES.TARGET.DandN.vcf"
			CMD=${CMD}" --databases_file_path ${ANNOVAR_DATABASE_FILE}"
			CMD=${CMD}" --ref_build_version ${ANNOVAR_REF_BUILD}"
			CMD=${CMD}" --threads ${THREADS}"
			CMD=${CMD}" --info_field_keys ${ANNOVAR_INFO_FIELD_KEYS_SPACED}"
			CMD=${CMD}" --header_mappings ${ANNOVAR_HEADER_MAPPINGS_SPACED}"
			CMD=${CMD}" --preserve_vcf_columns ${ANNOVAR_VCF_COLUMNS_SPACED}"
		CMD=${CMD}" --output_directory_path ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_ANNOVAR_TARGET/"
		CMD=${CMD}" &&"
		CMD=${CMD}" mv -v "
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_ANNOVAR_TARGET/${SM_TAG}.VARIANT_SITES.TARGET.DandN_ANNOVAR_REPORT.txt" \
			CMD=${CMD}" ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/ANNOVAR/${SM_TAG}.VARIANT_SITES.ON_TARGET_ANNOVAR_REPORT.txt"

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

END_ANNOVAR=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},T01,ANNOVAR,${HOSTNAME},${START_ANNOVAR},${END_ANNOVAR} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
