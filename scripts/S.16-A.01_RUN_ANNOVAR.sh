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

JAVA_1_6=$1
CIDRSEQSUITE_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6

# Copy the vcf that needs to have annovar run on it to the sample specific temp folder.

cp -rvf $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_TARGET/$SM_TAG".VARIANT_SITES.ON_TARGET.vcf" \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"

# Run annovar like a boss

START_RUN_ANNOVRAR=`date '+%s'`

$JAVA_1_6/java -jar $CIDRSEQSUITE_DIR/CIDRSeqSuite.jar \
-pipeline \
-annovar_directory_annotation \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR" \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"

END_RUN_ANNOVAR=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",Q.001,RUN_ANNOVAR,"$HOSTNAME","$START_RUN_ANNOVAR","$END_RUN_ANNOVAR \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_6/java -jar $CIDRSEQSUITE_DIR/CIDRSeqSuite.jar \
-pipeline \
-annovar_directory_annotation \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR" \
$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR" \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt
