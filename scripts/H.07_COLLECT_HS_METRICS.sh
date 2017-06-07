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
PICARD_DIR=$2
CORE_PATH=$3
SAMTOOLS_DIR=$4

PROJECT=$5
FAMILY=$6
SM_TAG=$7
REF_GENOME=$8
# BAIT_BED=$9
# TARGET_BED=${10}

# BAIT_NAME=`basename $BAIT_BED .bed`

# Create Picard style Calculate bed files (1-based start)

($SAMTOOLS_DIR/samtools view -H $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
| grep "@SQ" ; sed 's/\r//g' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_BAIT.bed" | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnBait.picard.bed"

($SAMTOOLS_DIR/samtools view -H $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
| grep "@SQ" ; sed 's/\r//g' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_TARGET.bed" | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTarget.picard.bed"

# NEED TO UPGRADE TO AN EVEN NEWER VERSION OF PICARD TO GET SOME OF THESE PARAMETERS...THAT I WANT

START_COLLECT_HS_METRICS=`date '+%s'`

$JAVA_1_8/java -jar $PICARD_DIR/picard.jar CollectHsMetrics \
INPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
PER_TARGET_COVERAGE=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE/$SM_TAG"_per_target_coverage.txt" \
REFERENCE_SEQUENCE=$REF_GENOME \
BAIT_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnBait.picard.bed" \
TARGET_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTarget.picard.bed" \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=10 \
BAIT_SET_NAME=Blah_bait_name \
VALIDATION_STRINGENCY=SILENT

END_COLLECT_HS_METRICS=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT"_BAM_REPORTS,Z.01,COLLECT_HS_METRICS,"$HOSTNAME","$START_COLLECT_HS_METRICS","$END_COLLECT_HS_METRICS \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $PICARD_DIR/picard.jar CollectHsMetrics \
INPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
PER_TARGET_COVERAGE=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE/$SM_TAG"_per_target_coverage.txt" \
REFERENCE_SEQUENCE=$REF_GENOME \
BAIT_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnBait.picard.bed" \
TARGET_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTarget.picard.bed" \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=10 \
BAIT_SET_NAME=blah_bait_name \
VALIDATION_STRINGENCY=SILENT \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE/$SM_TAG"_per_target_coverage.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
