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

## -----Haplotype Caller-----

## Call on Bait (padded or superset)

START_VARIANT_ANNOTATOR_GATHER=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz"

END_VARIANT_ANNOTATOR_GATHER=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",P.01-A.01,VARIANT_ANNOTATOR_GATHER,"$HOSTNAME","$START_VARIANT_ANNOTATOR_GATHER","$END_VARIANT_ANNOTATOR_GATHER \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".VQSR.ANNOTATED.vcf.gz.tbi" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
