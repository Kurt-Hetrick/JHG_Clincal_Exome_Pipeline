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
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
REF_GENOME=$7
KNOWN_INDEL_1=$8
KNOWN_INDEL_2=$9

## --Local realignment turn off downsampling

START_INDEL_REALIGNER=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
-R $REF_GENOME \
-known $KNOWN_INDEL_1 \
-known $KNOWN_INDEL_2 \
-targetIntervals $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/LOCAL_REALIGNMENT_INTERVALS/$SM_TAG"_LOCAL_REALIGNMENT_INTERVALS.intervals" \
-dt NONE \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".realign.bam"

END_INDEL_REALIGNER=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",E.01,INDEL_REALIGNER,"$HOSTNAME","$START_INDEL_REALIGNER","$END_INDEL_REALIGNER \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
-R $REF_GENOME \
-known $KNOWN_INDEL_1 \
-known $KNOWN_INDEL_2 \
-targetIntervals $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/LOCAL_REALIGNMENT_INTERVALS/$SM_TAG"_LOCAL_REALIGNMENT_INTERVALS.intervals" \
-dt NONE \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".realign.bam" \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"
