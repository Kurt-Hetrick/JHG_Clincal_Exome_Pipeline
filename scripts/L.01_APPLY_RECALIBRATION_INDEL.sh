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

START_APPLY_RECALIBRATION_INDEL=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.SNP.vcf" \
--ts_filter_level 99.9 \
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.tranches" \
-mode INDEL \
-o $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.SNP_INDEL.vcf"

END_APPLY_RECALIBRATION_INDEL=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",L.01,APPLY_RECALIBRATION_INDEL,"$HOSTNAME","$START_APPLY_RECALIBRATION_INDEL","$END_APPLY_RECALIBRATION_INDEL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.SNP.vcf" \
--ts_filter_level 99.9 \
-recalFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.recal" \
-tranchesFile $CORE_PATH/$PROJECT/$FAMILY/VCF/CONTROLS_PLUS_$FAMILY".HC.INDEL.tranches" \
-mode INDEL \
-o $CORE_PATH/$PROJECT/TEMP/CONTROLS_PLUS_$FAMILY".VQSR.SNP_INDEL.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
