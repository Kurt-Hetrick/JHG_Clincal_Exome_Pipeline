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
REF_GENOME=$6

# FILTER TO PASSING BIALLELIC SNVS IN THE CONTROLS PLUS FAMILY SET

START_FILTER_COHORT_SNV_PASS=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--excludeNonVariants \
--excludeFiltered \
-selectType SNP \
--restrictAllelesTo BIALLELIC \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz" \
-o $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.vcf"

END_FILTER_COHORT_SNV_PASS=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",S.01,FILTER_COHORT_SNV_ONLY_PASS_BIALLELEIC,"$HOSTNAME","$START_FILTER_COHORT_SNV_PASS","$END_FILTER_COHORT_SNV_PASS \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--excludeNonVariants \
--excludeFiltered \
-selectType SNP \
--restrictAllelesTo BIALLELIC \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz" \
-o $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
