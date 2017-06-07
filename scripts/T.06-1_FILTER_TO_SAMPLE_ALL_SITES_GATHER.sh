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

START_GATHER_SAMPLE_ALL_SITES=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.vcf.gz"

END_GATHER_SAMPLE_ALL_SITES=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",T.01,GATHER_SAMPLE_ALL_SITES_VCF,"$HOSTNAME","$START_GATHER_SAMPLE_ALL_SITES","$END_GATHER_SAMPLE_ALL_SITES \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.1.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.2.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.3.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.4.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.5.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.6.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.7.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.8.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.9.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.10.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.11.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.12.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.13.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.14.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.15.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.16.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.17.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.18.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.19.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.20.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.21.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.22.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.X.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".ALL_SITES.Y.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.vcf.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/FILTERED_ON_BAIT/$SM_TAG".ALL_SITES.vcf.gz.tbi" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
