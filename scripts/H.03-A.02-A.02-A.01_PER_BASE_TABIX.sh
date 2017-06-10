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

CORE_PATH=$1
TABIX_DIR=$2

PROJECT=$3
FAMILY=$4
SM_TAG=$5

START_PER_BASE_TABIX=`date '+%s'`

# Use TABIX to index padded per base refseq report

# index the per base report
# header is prefixed with #
# sequence name is column 1
# start and end is column 2

$TABIX_DIR/tabix -c \# -s 1 -b 2 -e 2 \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt.gz"

END_PER_BASE_TABIX=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,REFSEQ_PER_BASE_TABIX,"$HOSTNAME","$START_PER_BASE_TABIX","$END_PER_BASE_TABIX \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# echo awk $TABIX_DIR/tabix -c # -s 1 -b 2 -e 2 \
# $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt.gz" \
# >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt.gz.tbi" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
