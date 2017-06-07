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

JAVA_1_8=$1
PICARD_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
INPUT_BAM_FILE_STRING=$7

INPUT=`echo $INPUT_BAM_FILE_STRING | sed 's/,/ /g'`

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

## --Merge and Sort Bam files--

START_MERGE_BAM=`date '+%s'`

$JAVA_1_8/java -jar $PICARD_DIR/picard.jar MergeSamFiles \
$INPUT \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".original.bam" \
VALIDATION_STRINGENCY=SILENT \
SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true

END_MERGE_BAM=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",B.01,MERGE_BAM,"$HOSTNAME","$START_MERGE_BAM","$END_MERGE_BAM \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $PICARD_DIR/MergeSamFiles.jar \
$INPUT \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".original.bam" \
VALIDATION_STRINGENCY=SILENT \
SORT_ORDER=coordinate \
USE_THREADING=true \
CREATE_INDEX=true \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"
