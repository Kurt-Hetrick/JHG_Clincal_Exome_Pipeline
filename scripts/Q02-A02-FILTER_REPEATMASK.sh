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
	REF_GENOME=$5
	UCSC_REPEATMASK=$6
	MDUST_REPEATMASK=$7
	SAMPLE_SHEET=$8
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$9

# FILTER OUT REPEATMASKED REGIONS IN PASSING BIALLELIC SNVS IN FULL CALL SET

START_FILTER_COHORT_SNV_PASS_REPEATMASK=`date '+%s'`

	# construct command line

		CMD="singularity exec $ALIGNMENT_CONTAINER java -jar"
			CMD=$CMD" /gatk/gatk.jar"
		CMD=$CMD" SelectVariants"
			CMD=$CMD" --reference $REF_GENOME"
			CMD=$CMD" --exclude-intervals $UCSC_REPEATMASK"
			CMD=$CMD" --exclude-intervals $MDUST_REPEATMASK"
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/VCF_PREP/CONTROLS_PLUS_${FAMILY}.VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.vcf"
			CMD=$CMD" --output $CORE_PATH/$PROJECT/TEMP/VCF_PREP/CONTROLS_PLUS_${FAMILY}.VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.REPEATMASK.vcf.gz"

	# write command line to file and execute the command line

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/${FAMILY}_command_lines.txt
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/${FAMILY}_command_lines.txt
		echo $CMD | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

	# if exit does not equal 0 then exit with whatever the exit signal is at the end.
	# also write to file that this job failed

		if [ "$SCRIPT_STATUS" -ne 0 ]
		 then
			echo $FAMILY $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
			>> $CORE_PATH/$PROJECT/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
			exit $SCRIPT_STATUS
		fi

END_FILTER_COHORT_SNV_PASS=`date '+%s'`

# write out timing metrics to file

	echo ${FAMILY}_${PROJECT},S.01,FILTER_COHORT_SNV_PASS_REPEATMASK,${HOSTNAME},${START_FILTER_COHORT_SNV_PASS_REPEATMASK},${END_FILTER_COHORT_SNV_PASS_REPEATMASK} \
	>> $CORE_PATH/$PROJECT/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit $SCRIPT_STATUS