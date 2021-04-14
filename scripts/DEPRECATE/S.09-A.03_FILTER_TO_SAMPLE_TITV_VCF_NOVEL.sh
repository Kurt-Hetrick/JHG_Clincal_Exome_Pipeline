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
TITV_BED=$8
DBSNP_129=$9

# Filter to just on all of the variants all

START_FILTER_TO_SAMPLE_TITV_VCF_NOVEL=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
-L $TITV_BED \
--excludeNonVariants \
--excludeFiltered \
--discordance $DBSNP_129 \
--variant $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_BAIT/$SM_TAG".SNV.ON_BAIT.PASS.vcf.gz" \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".SNV.TITV_NOVEL.vcf"

END_FILTER_TO_SAMPLE_TITV_VCF_NOVEL=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",S.01,FILTER_TO_SAMPLE_TITV_VCF_NOVEL,"$HOSTNAME","$START_FILTER_TO_SAMPLE_TITV_VCF_NOVEL","$END_FILTER_TO_SAMPLE_TITV_VCF_NOVEL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
-L $TITV_BED \
--excludeNonVariants \
--excludeFiltered \
--discordance $DBSNP_129 \
--variant $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_BAIT/$SM_TAG".SNV.ON_BAIT.PASS.vcf.gz" \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".SNV.TITV_NOVEL.vcf" \
echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
