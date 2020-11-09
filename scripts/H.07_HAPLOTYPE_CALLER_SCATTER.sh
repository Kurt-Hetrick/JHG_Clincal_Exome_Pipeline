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
	SM_TAG=$4
	REF_GENOME=$5
	CODING_BED=$6
		CODING_BED_NAME=$(basename $CODING_BED .bed)
		CODING_MD5=$(md5sum $CODING_BED | cut -c 1-7)
	BAIT_BED=$7
		BAIT_BED_NAME=$(basename $BAIT_BED .bed)
	CHROMOSOME=$8
	GVCF_PAD=$9
	SAMPLE_SHEET=${10}
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=${11}

## -----Haplotype Caller BY CHROMOSOME -----
## CALL ON SUPER PADDED FILE
## CODING IS CONCATENATED WITH BAIT AND
## THEN PADDED WITH 250 BP AND THEN MERGED FOR OVERLAPPING INTERVALS

START_HAPLOTYPE_CALLER=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# construct command line

	CMD="singularity exec $GATK_3_7_0_CONTAINER java -jar" \
		CMD=$CMD" /usr/GenomeAnalysisTK.jar" \
	CMD=$CMD" -T HaplotypeCaller" \
		CMD=$CMD" -R $REF_GENOME" \
		CMD=$CMD" --input_file $CORE_PATH/$PROJECT/TEMP/$SM_TAG".bam"" \
		CMD=$CMD" -L $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME"-"$CODING_BED_NAME"-"$CODING_MD5"-"$GVCF_PAD"-BP-PAD-GVCF.bed"" \
		CMD=$CMD" -L $CHROMOSOME" \
		CMD=$CMD" --interval_set_rule INTERSECTION" \
		CMD=$CMD" --variant_index_type LINEAR" \
		CMD=$CMD" --variant_index_parameter 128000" \
		CMD=$CMD" -pairHMM VECTOR_LOGLESS_CACHING" \
		CMD=$CMD" --max_alternate_alleles 3" \
		CMD=$CMD" --emitRefConfidence BP_RESOLUTION" \
		CMD=$CMD" --annotation AS_BaseQualityRankSumTest" \
		CMD=$CMD" --annotation AS_FisherStrand" \
		CMD=$CMD" --annotation AS_InbreedingCoeff" \
		CMD=$CMD" --annotation AS_MappingQualityRankSumTest" \
		CMD=$CMD" --annotation AS_RMSMappingQuality" \
		CMD=$CMD" --annotation AS_ReadPosRankSumTest" \
		CMD=$CMD" --annotation AS_StrandOddsRatio" \
		CMD=$CMD" --annotation FractionInformativeReads" \
		CMD=$CMD" --annotation StrandBiasBySample" \
		CMD=$CMD" --annotation StrandAlleleCountsBySample" \
		CMD=$CMD" --annotation GCContent" \
		CMD=$CMD" --annotation AlleleBalanceBySample" \
		CMD=$CMD" --annotation AlleleBalance" \
		CMD=$CMD" --annotation LikelihoodRankSumTest" \
		CMD=$CMD" -bamout $CORE_PATH/$PROJECT/TEMP/$SM_TAG".HC."$CHROMOSOME".bam"" \
		CMD=$CMD" --emitDroppedReads" \
		CMD=$CMD" -o $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$CHROMOSOME".g.vcf.gz""

END_HAPLOTYPE_CALLER=`date '+%s'` # capture time process stops for wall clock tracking purposes.

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

# write out timing metrics to file

	echo $SM_TAG"_"$PROJECT",H.001,HAPLOTYPE_CALLER,"$HOSTNAME","$START_HAPLOTYPE_CALLER","$END_HAPLOTYPE_CALLER \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
