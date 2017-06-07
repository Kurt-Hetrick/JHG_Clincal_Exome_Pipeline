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

# Filter to just on all of the variants all

START_GATHER_FAMILY_ALL_SITES=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_TARGET/$SM_TAG".ALL_SITES.ON_TARGET.vcf.gz"

END_GATHER_FAMILY_ALL_SITES=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",T.01,GATHER_SAMPLE_ALL_SITES_TARGET_VCF,"$HOSTNAME","$START_GATHER_FAMILY_ALL_SITES","$END_GATHER_FAMILY_ALL_SITES \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.ON_TARGET.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_TARGET/$SM_TAG".ALL_SITES.ON_TARGET.vcf.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_TARGET/$SM_TAG".ALL_SITES.ON_TARGET.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_TARGET/$SM_TAG".ALL_SITES.ON_TARGET.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
