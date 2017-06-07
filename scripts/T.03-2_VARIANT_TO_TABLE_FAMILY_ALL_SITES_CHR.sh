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
CHROMOSOME=$7

# Filter to just on all of the variants all

START_VARIANT_TO_TABLE_FAMILY=`date '+%s'`

# not doing --splitMultiallelic here...maybe do one as an example and discuss with Molly
# do an example of molten output to look at/show molly

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantsToTable \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/TEMP/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY."$CHROMOSOME".vcf" \
--fields CHROM \
--fields POS \
--fields ID \
--fields REF \
--fields ALT \
--fields QUAL \
--fields EVENTLENGTH \
--fields FILTER \
--fields AC \
--fields AC_Orig \
--fields AN \
--fields AN_Orig \
--fields AF \
--fields AF_Orig \
--fields ABHet \
--fields ABHom \
--fields OND \
--fields HOM-REF \
--fields HET \
--fields HOM-VAR \
--fields NO-CALL \
--fields NSAMPLES \
--fields GQ_MEAN \
--fields GQ_STDDEV \
--fields NDA \
--fields Samples \
--fields DP \
--fields DP_Orig \
--fields QD \
--fields FS \
--fields SOR \
--fields MQ \
--fields MQRankSum \
--fields ReadPosRankSum \
--fields culprit \
--fields VQSLOD \
--fields POSITIVE_TRAIN_SITE \
--fields NEGATIVE_TRAIN_SITE \
--fields GC \
--fields FractionInformativeReads \
--fields HRun \
--fields RPA \
--fields RU \
--fields STR \
--fields MVLR \
--fields ClippingRankSum \
--fields OneKGP.AF \
--fields OneKGP.EAS_AF \
--fields OneKGP.AMR_AF \
--fields OneKGP.AFR_AF \
--fields OneKGP.EUR_AF \
--fields OneKGP.SAS_AF \
--genotypeFields AD \
--genotypeFields DP \
--genotypeFields GQ \
--genotypeFields GT \
--genotypeFields PGT \
--genotypeFields PID \
--genotypeFields PL \
--genotypeFields RGQ \
--genotypeFields SAC \
--allowMissingData \
--showFiltered \
-o $CORE_PATH/$PROJECT/TEMP/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY."$CHROMOSOME".txt"

END_VARIANT_TO_TABLE_FAMILY=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",Q.01,VARIANT_TO_TABLE_"$FAMILY"_ALL_SITES_"$CHROMOSOME","$HOSTNAME","$START_VARIANT_TO_TABLE_FAMILY","$END_VARIANT_TO_TABLE_FAMILY \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantsToTable \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/TEMP/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY."$CHROMOSOME".vcf" \
--fields CHROM \
--fields POS \
--fields ID \
--fields REF \
--fields ALT \
--fields QUAL \
--fields EVENTLENGTH \
--fields FILTER \
--fields AC \
--fields AC_Orig \
--fields AN \
--fields AN_Orig \
--fields AF \
--fields AF_Orig \
--fields ABHet \
--fields ABHom \
--fields OND \
--fields HOM-REF \
--fields HET \
--fields HOM-VAR \
--fields NO-CALL \
--fields NSAMPLES \
--fields GQ_MEAN \
--fields GQ_STDDEV \
--fields NDA \
--fields Samples \
--fields DP \
--fields DP_Orig \
--fields QD \
--fields FS \
--fields SOR \
--fields MQ \
--fields MQRankSum \
--fields ReadPosRankSum \
--fields culprit \
--fields VQSLOD \
--fields POSITIVE_TRAIN_SITE \
--fields NEGATIVE_TRAIN_SITE \
--fields GC \
--fields FractionInformativeReads \
--fields HRun \
--fields RPA \
--fields RU \
--fields STR \
--fields MVLR \
--fields ClippingRankSum \
--fields OneKGP.AF \
--fields OneKGP.EAS_AF \
--fields OneKGP.AMR_AF \
--fields OneKGP.AFR_AF \
--fields OneKGP.EUR_AF \
--fields OneKGP.SAS_AF \
--genotypeFields AD \
--genotypeFields DP \
--genotypeFields GQ \
--genotypeFields GT \
--genotypeFields PGT \
--genotypeFields PID \
--genotypeFields PL \
--genotypeFields RGQ \
--genotypeFields SAC \
--allowMissingData \
--showFiltered \
-o $CORE_PATH/$PROJECT/TEMP/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY."$CHROMOSOME".txt" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
