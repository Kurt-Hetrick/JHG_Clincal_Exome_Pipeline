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

# This would be a good candidate to write a bright module to load this.
source /isilon/cgc/PROGRAMS/change_R_version

set

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
REF_GENOME=$6
MILLS_1KG_GOLD_INDEL=$7

START_VARIANT_RECALIBRATOR_INDEL=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $REF_GENOME \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf" \
-resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_1KG_GOLD_INDEL \
--maxGaussians 4 \
--disable_auto_index_creation_and_locking_when_reading_rods \
-an QD \
-an FS \
-an SOR \
-an ReadPosRankSum \
-an MQRankSum \
-mode INDEL \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.8 \
-tranche 99.7 \
-tranche 99.6 \
-tranche 99.5 \
-tranche 99.4 \
-tranche 99.3 \
-tranche 99.2 \
-tranche 99.1 \
-tranche 99.0 \
-tranche 98.0 \
-tranche 97.0 \
-tranche 96.0 \
-tranche 95.0 \
-tranche 90.0 \
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.tranches" \
-rscriptFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.R"

END_VARIANT_RECALIBRATOR_INDEL=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",J.01,VARIANT_RECALIBRATOR_INDEL,"$HOSTNAME","$START_VARIANT_RECALIBRATOR_INDEL","$END_VARIANT_RECALIBRATOR_INDEL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $REF_GENOME \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf" \
-resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_1KG_GOLD_INDEL \
--maxGaussians 4 \
--disable_auto_index_creation_and_locking_when_reading_rods \
-an QD \
-an FS \
-an SOR \
-an ReadPosRankSum \
-an MQRankSum \
-mode INDEL \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.8 \
-tranche 99.7 \
-tranche 99.6 \
-tranche 99.5 \
-tranche 99.4 \
-tranche 99.3 \
-tranche 99.2 \
-tranche 99.1 \
-tranche 99.0 \
-tranche 98.0 \
-tranche 97.0 \
-tranche 96.0 \
-tranche 95.0 \
-tranche 90.0 \
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.tranches" \
-rscriptFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.R" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal.idx" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.tranches" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.R" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.R.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
