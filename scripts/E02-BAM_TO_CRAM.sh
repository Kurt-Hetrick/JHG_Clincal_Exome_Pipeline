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

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	REF_GENOME=$6
	THREADS=$7

	SAMPLE_SHEET=$8
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=$9

## --write lossless cram file. this is the deliverable

START_CRAM=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec ${ALIGNMENT_CONTAINER} samtools" \
		CMD=${CMD}" view" \
			CMD=${CMD}" -C ${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}.bam" \
			CMD=${CMD}" -T ${REF_GENOME}" \
			CMD=${CMD}" -@ ${THREADS}" \
			CMD=${CMD}" --write-index" \
			CMD=${CMD}" -O CRAM" \
		CMD=${CMD}" -o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram" \

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

END_CRAM=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},E01,CRAM,${HOSTNAME},${START_CRAM},${END_CRAM} \
	>> ${CORE_PATH}/${PROJECT}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

###############################################################################################################
## grab the read group header from the cram file and a make stub for the meta data ############################
## this is so if you need to regenerate a qc report and the cram file has been moved and not symlinked back ###
## this can be reference to use in generating the qc report. ##################################################
###############################################################################################################
#### THIS IS THE HEADER ################################################
#### "SM_TAG","PROJECT","RG_PU","LIBRARY" ##############################
#### "LIBRARY_PLATE","LIBRARY_WELL","LIBRARY_ROW","LIBRARY_COLUMN" #####
#### "HYB_PLATE","HYB_WELL","HYB_ROW","HYB_COLUMN" #####################
########################################################################

	if [ -f ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram ];
		then

			# grab field number for SM_TAG

				SM_FIELD=(`singularity exec ${ALIGNMENT_CONTAINER} samtools \
					view -H \
				${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^SM:/ \
						{print $1}'`)

			# grab field number for PLATFORM_UNIT_TAG

				PU_FIELD=(`singularity exec ${ALIGNMENT_CONTAINER} samtools \
					view -H \
				${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^PU:/ \
						{print $1}'`)

			# grab field number for LIBRARY_TAG

				LB_FIELD=(`singularity exec ${ALIGNMENT_CONTAINER} samtools \
					view -H \
				${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^LB:/ \
						{print $1}'`)

			# grab field number for PROGRAM_TAG

				PG_FIELD=(`singularity exec ${ALIGNMENT_CONTAINER} samtools \
					view -H \
				${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^PG:/ \
						{print $1}'`)

			# Now grab the header and format
				# breaking out the library name into its parts is assuming that the format is...
				# fill in empty fields with NA thing (for loop in awk) is a lifesaver
				# https://unix.stackexchange.com/questions/53448/replacing-missing-value-blank-space-with-zero

				singularity exec ${ALIGNMENT_CONTAINER} samtools \
					view -H \
				${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/CRAM/${SM_TAG}.cram \
					| grep ^@RG \
					| awk \
						-v SM_FIELD="$SM_FIELD" \
						-v PU_FIELD="$PU_FIELD" \
						-v LB_FIELD="$LB_FIELD" \
						-v PG_FIELD="$PG_FIELD" \
						'BEGIN {OFS="\t"} \
						{split($SM_FIELD,SMtag,":"); split($PU_FIELD,PU,":"); \
						split($LB_FIELD,Library,":"); split($PG_FIELD,Pipeline,":"); \
						print "'${PROJECT}'",SMtag[2],PU[2],Library[2],Pipeline[2]}' \
					| awk 'BEGIN { FS = OFS = "\t" } \
						{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' \
				>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt
		else
			echo -e "${PROJECT}\t${SM_TAG}\tNA\tNA\tNA" \
			| $DATAMASH_DIR/datamash transpose \
			>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt
	fi

# exit with the signal from samtools bam to cram

	exit ${SCRIPT_STATUS}
