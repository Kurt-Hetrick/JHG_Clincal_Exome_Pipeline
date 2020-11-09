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
	SM_TAG=$5
	REF_GENOME=$6
	BAIT_BED=$7
		BAIT_BED_NAME=(`basename $BAIT_BED .bed`)
	SAMPLE_SHEET=$8
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$9

## -----CONCATENATE SCATTERED g.vcf FILES INTO A SINGLE GRCh37 reference sorted g.vcf file-----

# Start with creating a *list file, reference sorted, to put into --variant.
# Assumption is that this is a correctly sorted GRCh37 reference file as the input reference used

	# Put the autosome into a file, sort numerically

		sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed" \
			| sed -r 's/[[:space:]]+/\t/g' \
			| cut -f 1 \
			| sort \
			| uniq \
			| awk '$1~/^[0-9]/' \
			| sort -k1,1n \
			| awk '{print "'$CORE_PATH'" "/" "'$PROJECT'" "/TEMP/" "'$SM_TAG'" "."$1".g.vcf.gz"}' \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gvcf.list"

	# Append X if present

		sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed" \
			| sed -r 's/[[:space:]]+/\t/g' \
			| cut -f 1 \
			| sort \
			| uniq \
			| awk '$1=="X"' \
			| awk '{print "'$CORE_PATH'" "/" "'$PROJECT'" "/TEMP/" "'$SM_TAG'" "."$1".g.vcf.gz"}' \
		>> $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gvcf.list"

	# Append Y if present

		sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed" \
			| sed -r 's/[[:space:]]+/\t/g' \
			| cut -f 1 \
			| sort \
			| uniq \
			| awk '$1=="Y"' \
			| awk '{print "'$CORE_PATH'" "/" "'$PROJECT'" "/TEMP/" "'$SM_TAG'" "."$1".g.vcf.gz"}' \
		>> $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gvcf.list"

## Gather the per chromosome gvcf files

START_HAPLOTYPE_CALLER_GATHER=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec $GATK_3_7_0_CONTAINER java -cp" \
			CMD=$CMD" /usr/GenomeAnalysisTK.jar" \
		CMD=$CMD" org.broadinstitute.gatk.tools.CatVariants" \
			CMD=$CMD" -R $REF_GENOME" \
			CMD=$CMD" --assumeSorted" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gvcf.list"" \
			CMD=$CMD" --outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.gz""

# write command line to file and execute the command line

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG"_command_lines.txt"
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

END_HAPLOTYPE_CALLER_GATHER=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo $SM_TAG"_"$PROJECT",H.01-A.01,HAPLOTYPE_CALLER_GATHER,"$HOSTNAME","$START_HAPLOTYPE_CALLER_GATHER","$END_HAPLOTYPE_CALLER_GATHER \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
