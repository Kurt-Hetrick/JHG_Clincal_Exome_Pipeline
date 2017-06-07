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

echo

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
REF_GENOME=$6
DBSNP=$7
HAPMAP=$8
OMNI_1KG=$9
HI_CONF_1KG_PHASE1_SNP=${10}

START_VARIANT_RECALIBRATOR_SNP=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $REF_GENOME \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf" \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
-resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1KG \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $HI_CONF_1KG_PHASE1_SNP \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
-mode SNP \
--disable_auto_index_creation_and_locking_when_reading_rods \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
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
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.tranches" \
-rscriptFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.R"

END_VARIANT_RECALIBRATOR_SNP=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",J.01,VARIANT_RECALIBRATOR_SNP,"$HOSTNAME","$START_VARIANT_RECALIBRATOR_SNP","$END_VARIANT_RECALIBRATOR_SNP \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $REF_GENOME \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf" \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
-resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1KG \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $HI_CONF_1KG_PHASE1_SNP \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
-mode SNP \
--disable_auto_index_creation_and_locking_when_reading_rods \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
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
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.tranches" \
-rscriptFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.R" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

###################
# Move the tranches PDF to TEMP so that it can be trashed.
mv $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.tranches.pdf" \
$CORE_PATH/$PROJECT/TEMP
###################

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.recal" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.recal.idx" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.tranches" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.R" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.SNV.R.pdf" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
