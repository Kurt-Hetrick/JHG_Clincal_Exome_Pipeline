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

CORE_PATH=$1
VERIFY_DIR=$2

PROJECT=$3
FAMILY=$4
SM_TAG=$5
CHROMOSOME=$6

## --Running verifyBamID--

START_VERIFYBAMID_CHR=`date '+%s'`

$VERIFY_DIR/verifyBamID \
--bam $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
--vcf $CORE_PATH/$PROJECT/TEMP/$SM_TAG".VerifyBamID."$CHROMOSOME".vcf" \
--out $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$CHROMOSOME \
--precise \
--verbose \
--maxDepth 2500

END_VERIFYBAMID_CHR=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT"_BAM_REPORTS,Z.09,VERIFYBAMID_"$CHROMOSOME","$HOSTNAME","$START_VERIFYBAMID_CHR","$END_VERIFYBAMID_CHR \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $VERIFY_DIR/verifyBamID \
--bam $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
--vcf $CORE_PATH/$PROJECT/TEMP/$SM_TAG".VerifyBamID."$CHROMOSOME".vcf" \
--out $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$CHROMOSOME \
--precise \
--verbose \
--maxDepth 2500 \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"
