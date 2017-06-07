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
GENE_LIST=$4

PROJECT=$5
FAMILY=$6
SM_TAG=$7
REF_GENOME=$8
# TARGET_BED=$9

### --Depth of Coverage On Target--

START_DOC_TARGET=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_TARGET.bed" \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED" \
-ct 10 \
-ct 15 \
-ct 20 \
-ct 30 \
-ct 50

END_DOC_TARGET=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT"_BAM_REPORTS,Z.01,DOC_TARGET,"$HOSTNAME","$START_DOC_TARGET","$END_DOC_TARGET \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_TARGET.bed" \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-omitBaseOutput \
-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED" \
-ct 10 \
-ct 15 \
-ct 20 \
-ct 30 \
-ct 50 \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.EveryBase.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.EveryBase.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_counts" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_counts.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_counts.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_proportions" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_proportions.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_cumulative_coverage_proportions.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_gene_summary" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_gene_summary.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_gene_summary.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_statistics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_statistics.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_statistics.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_summary" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_summary.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_interval_summary.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_statistics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_statistics.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_statistics.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_summary" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_summary.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/TARGET/$SM_TAG".TARGET.BED.sample_summary.csv" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
