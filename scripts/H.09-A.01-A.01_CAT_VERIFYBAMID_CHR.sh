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
CORE_PATH=$1

PROJECT=$2
FAMILY=$3
SM_TAG=$4

echo \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG.verifybamID_joined.txt

for CHROMOSOME in {1..22}
do
cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$CHROMOSOME".selfSM" \
| grep -v ^# \
| awk 'BEGIN {OFS="\t"} {print($1,"'$CHROMOSOME'",$7,$4,$8,$9,$6)}' \
>> $CORE_PATH/$PROJECT/TEMP/$SM_TAG.verifybamID_joined.txt
done

sed -i '/^\s*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG.verifybamID_joined.txt

echo "#SM_TAG" CHROM VERIFYBAM_FREEMIX VERIFYBAM_SNPS VERIFYBAM_FREELK1 VERRIFYBAM_FREELK0 VERIFYBAM_AVG_DP \
>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID_CHR/$SM_TAG.VERIFYBAMID.PER_CHR.txt

cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG.verifybamID_joined.txt \
>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID_CHR/$SM_TAG.VERIFYBAMID.PER_CHR.txt

sed -i 's/ /\t/g' $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID_CHR/$SM_TAG.VERIFYBAMID.PER_CHR.txt

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID_CHR/$SM_TAG.VERIFYBAMID.PER_CHR.txt \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
