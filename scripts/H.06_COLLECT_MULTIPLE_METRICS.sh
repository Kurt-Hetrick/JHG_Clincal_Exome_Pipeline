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
DBSNP=$9
TI_TV_BED=${10}

# blah, blah, blah

# Create a picard style target bed file. This is used for CollectSequencingArtifactMetrics...it should not be used for anything else...i think

($SAMTOOLS_DIR/samtools view -H $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
| grep "@SQ" ; sed 's/\r//g' $TI_TV_BED | awk '{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' | sed 's/ /\t/g') \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTiTv.picard_2.bed"

START_COLLECT_MULTIPLE_METRICS=`date '+%s'`

$JAVA_1_8/java -jar $PICARD_DIR/picard.jar \
CollectMultipleMetrics \
INPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$SM_TAG \
REFERENCE_SEQUENCE=$REF_GENOME \
DB_SNP=$DBSNP \
INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTiTv.picard_2.bed" \
PROGRAM=CollectGcBiasMetrics \
PROGRAM=CollectSequencingArtifactMetrics

END_COLLECT_MULTIPLE_METRICS=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT"_BAM_REPORTS,Z.01,COLLECT_MULTIPLE_METRICS,"$HOSTNAME","$START_COLLECT_MULTIPLE_METRICS","$END_COLLECT_MULTIPLE_METRICS \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $PICARD_DIR/picard.jar \
CollectMultipleMetrics \
INPUT=$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$SM_TAG \
REFERENCE_SEQUENCE=$REF_GENOME \
DB_SNP=$DBSNP \
INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTiTv.picard_2.bed" \
PROGRAM=CollectGcBiasMetrics \
PROGRAM=CollectSequencingArtifactMetrics \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

# Move and rename bait bais metrics/summary files to the reports directory and add a txt extension

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".bait_bias_detail_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/METRICS/$SM_TAG".bait_bias_detail_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/METRICS/$SM_TAG".bait_bias_detail_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".bait_bias_summary_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# Move and rename pre adapter metrics/summary files to the reports directory and add a txt extension

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".pre_adapter_detail_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/METRICS/$SM_TAG".pre_adapter_detail_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/METRICS/$SM_TAG".pre_adapter_detail_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".pre_adapter_summary_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the results from collect alignment summary metrics to the reports folder

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".alignment_summary_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the base distribution by cycle reports into the report directory

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".base_distribution_by_cycle.pdf" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/PDF/$SM_TAG".base_distribution_by_cycle.pdf"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/PDF/$SM_TAG".base_distribution_by_cycle.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".base_distribution_by_cycle_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the insert size reports into the report directory

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".insert_size_histogram.pdf" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/PDF/$SM_TAG".insert_size_histogram.pdf"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/PDF/$SM_TAG".insert_size_histogram.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".insert_size_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the mean quality by cycle into the report directory

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".quality_by_cycle.pdf" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/PDF/$SM_TAG".quality_by_cycle.pdf"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/PDF/$SM_TAG".quality_by_cycle.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".quality_by_cycle_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/METRICS/$SM_TAG".quality_by_cycle_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/METRICS/$SM_TAG".quality_by_cycle_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the basecall by q score into the report directory

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".quality_distribution.pdf" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/PDF/$SM_TAG".quality_distribution.pdf"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/PDF/$SM_TAG".quality_distribution.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".quality_distribution_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/METRICS/$SM_TAG".quality_distribution_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/METRICS/$SM_TAG".quality_distribution_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

# move the gc bias reports into the report directory

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gc_bias.pdf" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/PDF/$SM_TAG".gc_bias.pdf"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/PDF/$SM_TAG".gc_bias.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gc_bias.detail_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/METRICS/$SM_TAG".gc_bias.detail_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/METRICS/$SM_TAG".gc_bias.detail_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

mv -v $CORE_PATH/$PROJECT/TEMP/$SM_TAG".gc_bias.summary_metrics" \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/SUMMARY/$SM_TAG".gc_bias.summary_metrics.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/SUMMARY/$SM_TAG".gc_bias.summary_metrics.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
