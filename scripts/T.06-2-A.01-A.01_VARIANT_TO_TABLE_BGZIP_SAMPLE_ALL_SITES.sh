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

TABIX_DIR=$1
CORE_PATH=$2

PROJECT=$3
FAMILY=$4
SM_TAG=$5

# Filter to just on all of the variants all

START_VARIANT_TO_TABLE_BGZIP_SAMPLE=`date '+%s'`

# not doing --splitMultiallelic here...maybe do one as an example and discuss with Molly
# do an example of molten output to look at/show molly


$TABIX_DIR/bgzip \
-c $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.txt" \
>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.txt.gz"

# rm -rvf /usr/local/sandbox/cgc/$SM_TAG".ALL_SITES.txt"

END_VARIANT_TO_TABLE_BGZIP_SAMPLE=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",U.01,VARIANT_TO_TABLE_BGZIP_"$SM_TAG"_ALL_SITES,"$HOSTNAME","$START_VARIANT_TO_TABLE_BGZIP_SAMPLE","$START_VARIANT_TO_TABLE_BGZIP_SAMPLE \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $TABIX_DIR/bgzip \
-c /usr/local/sandbox/cgc/$SM_TAG".ALL_SITES.txt" \
\>\| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.txt.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.txt.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
