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

START_GENOTYPE_GVCF_GATHER=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf"

END_GENOTYPE_GVCF_GATHER=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",I.01-A.01,GENOTYPE_GVCF_GATHER,"$HOSTNAME","$START_GENOTYPE_GVCF_GATHER","$END_GENOTYPE_GVCF_GATHER \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".RAW.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
