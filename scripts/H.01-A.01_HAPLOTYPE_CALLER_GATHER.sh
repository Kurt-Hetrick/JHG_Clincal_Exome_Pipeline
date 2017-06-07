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

## -----Haplotype Caller-----

## Call on Bait (padded or superset)

START_HAPLOTYPE_CALLER_GATHER=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".1.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".2.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".3.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".4.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".5.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".6.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".7.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".8.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".9.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".10.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".11.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".12.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".13.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".14.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".15.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".16.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".17.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".18.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".19.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".20.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".21.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".22.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".X.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".Y.g.vcf" \
--outputFile $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.gz"

END_HAPLOTYPE_CALLER_GATHER=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.01-A.01,HAPLOTYPE_CALLER_GATHER,"$HOSTNAME","$START_HAPLOTYPE_CALLER_GATHER","$END_HAPLOTYPE_CALLER_GATHER \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".1.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".2.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".3.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".4.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".5.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".6.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".7.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".8.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".9.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".10.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".11.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".12.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".13.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".14.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".15.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".16.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".17.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".18.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".19.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".20.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".21.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".22.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".X.g.vcf" \
--variant $CORE_PATH/$PROJECT/TEMP/$SM_TAG".Y.g.vcf" \
-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.gz" \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF/$SM_TAG".g.vcf.gz.tbi" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
