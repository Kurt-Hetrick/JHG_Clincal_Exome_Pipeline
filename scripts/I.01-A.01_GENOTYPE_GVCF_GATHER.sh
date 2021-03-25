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
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$7

# GATHER UP PER CHROMOSOME GENOTYPED VCF FILES.

START_GENOTYPE_GVCF_GATHER=`date '+%s'`

	# construct command line

		CMD="singularity exec $GATK_3_7_0_CONTAINER java -cp" \
				CMD=$CMD" /usr/GenomeAnalysisTK.jar" \
		CMD=$CMD" org.broadinstitute.gatk.tools.CatVariants" \
			CMD=$CMD" -R $REF_GENOME" \
			CMD=$CMD" --assumeSorted" \

			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.1.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.2.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.3.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.4.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.5.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.6.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.7.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.8.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.9.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.10.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.11.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.12.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.13.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.14.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.15.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.16.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.17.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.18.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.19.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.20.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.21.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.22.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.X.vcf"" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.Y.vcf"" \

			CMD=$CMD" --outputFile $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf""

	# write command line to file and execute the command line

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/$FAMILY"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$FAMILY"_command_lines.txt"
		echo $CMD | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

	# if exit does not equal 0 then exit with whatever the exit signal is at the end.
	# also write to file that this job failed

		if [ "$SCRIPT_STATUS" -ne 0 ]
		 then
			echo $SM_TAG $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
			>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
			exit $SCRIPT_STATUS
		fi

END_GENOTYPE_GVCF_GATHER=`date '+%s'`

# write out timing metrics to file

	echo $FAMILY"_"$PROJECT",I.01-A.01,GENOTYPE_GVCF_GATHER,"$HOSTNAME","$START_GENOTYPE_GVCF_GATHER","$END_GENOTYPE_GVCF_GATHER \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
