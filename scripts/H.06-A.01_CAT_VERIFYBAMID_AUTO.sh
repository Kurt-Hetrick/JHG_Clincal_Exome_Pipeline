# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

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

# INPUT VARIABLES

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	BAIT_BED=$6
		BAIT_BED_NAME=(`basename $BAIT_BED .bed`)

# CREATE A BLANK FILE TO START ADDING TO

	echo \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_noheader.txt"

# LOOP THROUGH THE PER AUTOSOME VERIFYBAMID OUTPUT AND ADD THEM TO ABOVE EMPTY FILE

	for AUTOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed" \
				| sed -r 's/[[:space:]]+/\t/g' \
				| cut -f 1 \
				| egrep -v "X|Y|MT" \
				| sort -k 1,1n \
				| uniq \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					collapse 1 \
				| sed 's/,/ /g');
		do
			# I'm stripping out the "chr" prefix here b/c I don't want to deal with it...I should specify the column in case SM_TAG contain chr...
			cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$AUTOSOME".selfSM" \
				| grep -v ^# \
				| awk 'BEGIN {OFS="\t"} {print($1,"'$AUTOSOME'",$7,$4,$8,$9,$6)}' \
				| sed 's/chr//g' \
			>> $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_noheader.txt"
	done

# REMOVE BLANK LINES

	sed -i '/^\s*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_noheader.txt"

# ADD HEADER

	echo "#SM_TAG" CHROM VERIFYBAM_FREEMIX VERIFYBAM_SNPS VERIFYBAM_FREELK1 VERRIFYBAM_FREELK0 VERIFYBAM_AVG_DP \
		| cat - $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_noheader.txt" \
	>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID_AUTO/$SM_TAG".VERIFYBAMID.PER_AUTOSOME.txt"
