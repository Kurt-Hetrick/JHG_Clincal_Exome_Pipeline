# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q cgc.q

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

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
REF_GENOME=$7
CHROMOSOME=$8

# Filter to just on all of the variants all

START_FILTER_TO_SAMPLE_ALL_SITES=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
--keepOriginalDP \
--removeUnusedAlternates \
--sample_name $SM_TAG \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED."$CHROMOSOME".vcf" \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES."$CHROMOSOME".vcf"

END_FILTER_TO_SAMPLE_ALL_SITES=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",S.01,FILTER_TO_SAMPLE_ALL_SITES"$CHROMOSOME","$HOSTNAME","$START_FILTER_TO_SAMPLE_ALL_SITES","$END_FILTER_TO_SAMPLE_ALL_SITES \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
--keepOriginalDP \
--removeUnusedAlternates \
--sample_name $SM_TAG \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED."$CHROMOSOME".vcf" \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES."$CHROMOSOME".vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
