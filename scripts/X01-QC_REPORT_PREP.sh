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

	set

	echo

# INPUT VARIABLES

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5
	FATHER=$6
	MOTHER=$7
	GENDER=$8
	PHENOTYPE=$9

#########################################################
##### Grabbing the BAM header (for RG ID,PU,LB,etc) #####
################################################################################
##### THIS IS THE HEADER #######################################################
##### "PROJECT","SM_TAG","PLATFORM_UNIT","LIBRARY_NAME","PIPELINE_VERSION" #####
##### "FAMILY","FATHER","MOTHER","LIMS_SEX","PHENOTYPE" ########################
################################################################################

	if [ -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt ]
		then
			cat $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					-s \
					-g 1,2 \
					collapse 3 \
					unique 4 \
					unique 5 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} \
					{print $0,"'$FAMILY'","'$FATHER'","'$MOTHER'","'$GENDER'","'$PHENOTYPE'"}' \
				| awk 'BEGIN {OFS="\t"} \
					$9=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,"MALE",$10} \
					$9=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"FEMALE",$10} \
					$9!="1"&&$9!="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"UNKNOWN",$10}' \
				| awk 'BEGIN {OFS="\t"} \
					$10=="-9" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \
					$10=="0" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \
					$10=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"UNAFFECTED"} \
					$10=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"AFFECTED"}' \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					transpose \
			>| $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		elif [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt && \
			-f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram ]];
			then

			# grab field number for SM_TAG

				SM_FIELD=(`singularity exec $ALIGNMENT_CONTAINER samtools view -H \
				$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^SM:/ {print $1}'`)

			# grab field number for PLATFORM_UNIT_TAG

				PU_FIELD=(`singularity exec $ALIGNMENT_CONTAINER samtools view -H \
				$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^PU:/ {print $1}'`)

			# grab field number for LIBRARY_TAG

				LB_FIELD=(`singularity exec $ALIGNMENT_CONTAINER samtools view -H \
				$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^LB:/ {print $1}'`)

			# grab field number for PROGRAM_TAG

				PG_FIELD=(`singularity exec $ALIGNMENT_CONTAINER samtools view -H \
				$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram \
					| grep -m 1 ^@RG \
					| sed 's/\t/\n/g' \
					| cat -n \
					| sed 's/^ *//g' \
					| awk '$2~/^PG:/ {print $1}'`)

			# Now grab the header and format
				# fill in empty fields with NA thing (for loop in awk) is a lifesaver
				# https://unix.stackexchange.com/questions/53448/replacing-missing-value-blank-space-with-zero

				singularity exec $ALIGNMENT_CONTAINER samtools \
					view -H \
				$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM/${SM_TAG}.cram \
					| grep ^@RG \
					| awk \
						-v SM_FIELD="$SM_FIELD" \
						-v PU_FIELD="$PU_FIELD" \
						-v LB_FIELD="$LB_FIELD" \
						-v PG_FIELD="$PG_FIELD" \
						'BEGIN {OFS="\t"} \
						{split($SM_FIELD,SMtag,":"); \
						split($PU_FIELD,PU,":"); \
						split($LB_FIELD,Library,":"); \
						split($PG_FIELD,Pipeline,":"); \
						print "'$PROJECT'",SMtag[2],PU[2],Library[2],Pipeline[2]}' \
					| awk 'BEGIN { FS = OFS = "\t" } \
						{ for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' \
					| singularity exec $ALIGNMENT_CONTAINER datamash \
						-s \
						-g 1,2 \
						collapse 3 \
						unique 4 \
						unique 5 \
					| sed 's/,/;/g' \
					| awk 'BEGIN {OFS="\t"} \
						{print $0,"'$FAMILY'","'$FATHER'","'$MOTHER'","'$GENDER'","'$PHENOTYPE'"}' \
					| awk 'BEGIN {OFS="\t"} \
						$9=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,"MALE",$10} \
						$9=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"FEMALE",$10} \
						$9!="1"&&$9!="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"UNKNOWN",$10}' \
					| awk 'BEGIN {OFS="\t"} \
						$10=="-9" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \
						$10=="0" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \
						$10=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"UNPHENOTYPE"} \
						$10=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"PHENOTYPE"}' \
					| singularity exec $ALIGNMENT_CONTAINER datamash \
						transpose \
				>| $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
		else
			echo -e "$PROJECT\t$SM_TAG\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					transpose \
			>| $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

#################################################
##### GENDER CHECK FROM ANEUPLOIDY CHECK ########
#################################################
##### THIS IS THE HEADER ########################
##### X_AVG_DP,X_NORM_DP,Y_AVG_DP,Y_NORM_DP #####
#################################################

	awk 'BEGIN {OFS="\t"} \
		$2=="X"&&$3=="whole" {print "X",$6,$7} \
		$2=="Y"&&$3=="whole" {print "Y",$6,$7}' \
	$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK/${SM_TAG}.chrom_count_report.txt \
		| paste - - \
		| awk 'BEGIN {OFS="\t"} \
			END {if ($1=="X"&&$4=="Y") print $2,$3,$5,$6 ; \
			else if ($1=="X"&&$4=="") print $2,$3,"NaN","NaN" ; \
			else if ($1=="Y"&&$4=="") print "NaN","NaN",$5,$6 ; \
			else print "NaN","NaN","NaN","NaN"}' \
		| singularity exec $ALIGNMENT_CONTAINER datamash \
			transpose \
	>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

#########################################################################################
##### VERIFY BAM ID #####################################################################
#########################################################################################
##### THIS IS THE HEADER ################################################################
##### "VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0" #####
##### "VERIFYBAM_DIFF_LK0_LK1","VERIFYBAM_AVG_DP" #######################################
#########################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID/${SM_TAG}.selfSM ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {OFS="\t"} \
				NR>1 \
				{print $7*100,$4,$8,$9,($9-$8),$6}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID/${SM_TAG}.selfSM \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

####################################################################################
##### INSERT SIZE ##################################################################
####################################################################################
##### THIS IS THE HEADER ###########################################################
##### "MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE" #####
####################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/${SM_TAG}.insert_size_metrics.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {OFS="\t"} \
				NR==8 \
				{print $1,$6,$7}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/${SM_TAG}.insert_size_metrics.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

##########################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 1 ###############################
##########################################################################
##### THIS THE HEADER ####################################################
##### "PCT_PF_READS_ALIGNED_R1","PF_HQ_ALIGNED_READS_R1" #################
##### "PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1" #####
##### "PCT_READS_ALIGNED_IN_PAIRS_R1","PCT_ADAPTER_R1" ###################
##########################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {OFS="\t"} \
				NR==8 \
				{if ($1=="UNPAIRED") print "0","0","0","0","0","0","0"; \
				else print $7*100,$9,$13,$14,$15,$18*100,$24*100}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

##########################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 2 ###############################
##########################################################################
##### THIS THE HEADER ####################################################
##### "PCT_PF_READS_ALIGNED_R2","PF_HQ_ALIGNED_READS_R2" #################
##### "PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2" #####
##### "PCT_READS_ALIGNED_IN_PAIRS_R2","PCT_ADAPTER_R2" ###################
##########################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {OFS="\t"} \
				NR==9 \
				{if ($1=="") print "0","0","0","0","0","0","0"; \
				else print $7*100,$9,$13,$14,$15,$18*100,$24*100}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

#######################################################################################
##### ALIGNMENT SUMMARY METRICS FOR PAIR ##############################################
#######################################################################################
##### THIS THE HEADER ####################################################################
##### "TOTAL_READS","RAW_GIGS","PCT_PF_READS_ALIGNED_PAIR","PF_MISMATCH_RATE_PAIR" #######
##### "PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR","PCT_READS_ALIGNED_IN_PAIRS_PAIR" #####
##### "PCT_PF_READS_IMPROPER_PAIRS_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR" #######
##########################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {OFS="\t"} \
				NR==10 \
				{if ($1=="") print "0","0","0","0","0","0","0","0","0","0" ; \
				else print $2,($2*$16/1000000000),$7*100,$13,$14,$15,$18*100,$20*100,$22,$23*100}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/${SM_TAG}.alignment_summary_metrics.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

####################################################################################
##### MARK DUPLICATES REPORT #######################################################
####################################################################################
##### THIS IS THE HEADER ###########################################################
##### "UNMAPPED_READS","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION" ########
##### "ESTIMATED_LIBRARY_SIZE","SECONDARY_OR_SUPPLEMENTARY_READS" ##################
##### "READ_PAIR_DUPLICATES","READ_PAIRS_EXAMINED","PAIRED_DUP_RATE" ###############
##### "UNPAIRED_READ_DUPLICATES","UNPAIRED_READS_EXAMINED","UNPAIRED_DUP_RATE" #####
##### "PERCENT_DUPLICATION_OPTICAL" ################################################
####################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PICARD_DUPLICATES/${SM_TAG}_MARK_DUPLICATES.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
		else

			MAX_RECORD=(`grep -n "^$" \
					$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PICARD_DUPLICATES/${SM_TAG}_MARK_DUPLICATES.txt \
					| awk 'BEGIN {FS=":"} NR==2 {print $1}'`)

			awk 'BEGIN {OFS="\t"} \
				NR>7&&NR<'$MAX_RECORD' \
				{if ($10!~/[0-9]/) print $5,$8,"NaN","NaN",$4,$7,$3,"NaN",$6,$2,"NaN" ; \
				else if ($10~/[0-9]/&&$2=="0") print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,"NaN" ; \
				else print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,($6/$2)}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PICARD_DUPLICATES/${SM_TAG}_MARK_DUPLICATES.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				sum 1 \
				sum 2 \
				mean 4 \
				sum 5 \
				sum 6 \
				sum 7 \
				sum 9 \
				sum 10 \
			| awk 'BEGIN {OFS="\t"} \
				{if ($3!~/[0-9]/) print $1,$2,"NaN","NaN",$4,$5,$6,"NaN",$7,$8,"NaN","NaN" ; \
				else if ($3~/[0-9]/&&$1=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else if ($3~/[0-9]/&&$1!="0"&&$8=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,($7/$8),($2/$6)*100}' \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

#######################################################################################################
##### HYBRIDIZATION SELECTION REPORT ##################################################################
#######################################################################################################
##### THIS IS THE HEADER ##############################################################################
##### "GENOME_SIZE","BAIT_SET","BAIT_TERRITORY","TARGET_TERRITORY" ####################################
##### "PCT_PF_UQ_READS_ALIGNED","PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","ON_BAIT_VS_SELECTED" #######
##### "MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE","MAX_TARGET_COVERAGE" ######
##### "PCT_USABLE_BASES_ON_BAIT","ZERO_CVG_TARGETS_PCT" ###############################################
##### "PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET" ###########################
##### "PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X" #####
##### "AT_DROPOUT","GC_DROPOUT","THEORETICAL_HET_SENSITIVITY","HET_SNP_Q" #############################
#######################################################################################################

	# this will take when there are no reads in the file...but i don't think that it will handle when there are reads, but none fall on target
	# the next time i that happens i'll fix this to handle it.

		if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/${SM_TAG}_hybridization_selection_metrics.txt ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			awk 'BEGIN {FS="\t";OFS="\t"} \
				NR==8 \
				{if ($12=="?"&&$44=="") \
					print $2,$1,$3,$4,"NaN",($14/1000000000),"NaN","NaN",$22,$23,$24,$25,"NaN",$29,"NaN","NaN","NaN","NaN",$39,$40,$41,$42,$51,$52,$53,$54 ; \
				else if ($12!="?"&&$44=="") \
					print $2,$1,$3,$4,$12*100,($14/1000000000),$19*100,$21,$22,$23,$24,$25,$26*100,$29*100,$31*100,$32*100,$33*100,$34*100,$39*100,$40*100,$41*100,$42*100,$51,$52,$53,$54 ; \
				else print $2,$1,$3,$4,$12*100,($14/1000000000),$19*100,$21,$22,$23,$24,$25,$26*100,$29*100,$31*100,$32*100,$33*100,$34*100,$39*100,$40*100,$41*100,$42*100,$51,$52,$53,$54}' \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/${SM_TAG}_hybridization_selection_metrics.txt \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
		fi

##############################################
##### BAIT BIAS REPORT FOR Cref and Gref #####
##############################################
##### THIS IS THE HEADER #####################
##### "Cref_Q","Gref_Q" ######################
##############################################

	if [[ ! -f $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/${SM_TAG}.bait_bias_summary_metrics.txt ]]
		then
			echo -e NaN'\t'NaN \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt

		else
			grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/${SM_TAG}.bait_bias_summary_metrics.txt \
				| sed '/^$/d' \
				| awk 'BEGIN {OFS="\t"} $12=="Cref"||$12=="Gref" {print $5}' \
				| paste - - \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					collapse 1 \
					collapse 2 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} {print $0}' \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					transpose \
			>> $CORE_PATH/$PROJECT/TEMP/${SM_TAG}.QC_REPORT_TEMP.txt
	fi

# ############################################################
# ##### PRE-ADAPTER BIAS REPORT FOR Deamination and OxoG #####
# ############################################################
# ##### THIS IS THE HEADER ###################################
# ##### SM_TAG,Deamination_Q,OxoG_Q ##########################
# ############################################################

# grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" \
# | sed '/^$/d' \
# | awk 'BEGIN {OFS="\t"} $12=="Deamination"||$12=="OxoG"  {print $5}' \
# | paste - - \
# | awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_PRE_ADAPTER.TXT"

# ###########################################################################
# ##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT SNVS ########################
# ###########################################################################
# ##### THIS IS THE HEADER ##################################################
# ##### "SM_TAG""\t""COUNT_SNV_ON_BAIT""\t""PERCENT_SNV_ON_BAIT_SNP138" ##### 
# ###########################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_BAIT/$SM_TAG".SNV.ON_BAIT.PASS.vcf.gz" \
# | awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
# END {if (SNV_COUNT!="") {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
# else {print "'$SM_TAG'","0","NaN"}}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_SNV_METRICS.TXT"

# ###############################################################################
# ##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET SNVS ##########################
# ###############################################################################
# ##### THIS IS THE HEADER ######################################################
# ##### "SM_TAG""\t""COUNT_SNV_ON_TARGET""\t""PERCENT_SNV_ON_TARGET_SNP138" ##### 
# ###############################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_TARGET/$SM_TAG".SNV.ON_TARGET.PASS.vcf.gz" \
# | awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
# END {if (SNV_COUNT!="") {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
# else {print "'$SM_TAG'","0","NaN"}}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_SNV_METRICS.TXT"

# ##############################################################
# ##### GRABBING TI/TV ON UCSC CODING EXONS, ALL ###############
# ##############################################################
# ##### THIS IS THE HEADER #####################################
# ##### "SM_TAG""\t""ALL_TI_TV_COUNT""\t""ALL_TI_TV_RATIO" #####
# ##############################################################

# awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
# else {print "'$SM_TAG'","0","NaN"}}' \
# $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_All_.titv.txt" \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_ALL.TXT"

# ##################################################################
# ##### GRABBING TI/TV ON UCSC CODING EXONS, KNOWN #################
# ##################################################################
# ##### THIS IS THE HEADER #########################################
# ##### "SM_TAG""\t""KNOWN_TI_TV_COUNT""\t""KNOWN_TI_TV_RATIO" #####
# ##################################################################

# awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
# else {print "'$SM_TAG'","0","NaN"}}' \
# $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Known_.titv.txt" \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_KNOWN.TXT"

# ##################################################################
# ##### GRABBING TI/TV ON UCSC CODING EXONS, NOVEL #################
# ##################################################################
# ##### THIS IS THE HEADER #########################################
# ##### "SM_TAG""\t""NOVEL_TI_TV_COUNT""\t""NOVEL_TI_TV_RATIO" #####
# ##################################################################

# awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
# else {print "'$SM_TAG'","0","NaN"}}' \
# $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt" \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_NOVEL.TXT"

# ######################################################################################################################################
# ##### INDEL METRICS ON BAIT ##########################################################################################################
# ######################################################################################################################################
# ##### THIS IS THE HEADER #############################################################################################################
# ##### "SM_TAG","COUNT_ALL_INDEL_BAIT","ALL_INDEL_BAIT_PCT_SNP138","COUNT_BIALLELIC_INDEL_BAIT","BIALLELIC_INDEL_BAIT_PCT_SNP138" #####
# ######################################################################################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_BAIT/$SM_TAG".INDEL.ON_BAIT.PASS.vcf.gz" \
# | awk '{INDEL_COUNT++NR} \
# {INDEL_BIALLELIC+=($5!~",")} \
# {DBSNP_COUNT+=($3~"rs")} \
# {DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
# END {if (INDEL_BIALLELIC==""&&INDEL_COUNT=="") print "'$SM_TAG'","0","NaN","0","NaN"; \
# else if (INDEL_BIALLELIC==0&&INDEL_COUNT>=1) print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN"; \
# else print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_INDEL_METRICS.TXT"

# ##############################################################################################################################################
# ##### INDEL METRICS ON TARGET ################################################################################################################
# ##############################################################################################################################################
# ##### THIS IS THE HEADER #####################################################################################################################
# ##### "SM_TAG","COUNT_ALL_INDEL_TARGET","ALL_INDEL_TARGET_PCT_SNP138","COUNT_BIALLELIC_INDEL_TARGET","BIALLELIC_INDEL_TARGET_PCT_SNP138" #####
# ##############################################################################################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_TARGET/$SM_TAG".INDEL.ON_TARGET.PASS.vcf.gz" \
# | awk '{INDEL_COUNT++NR} \
# {INDEL_BIALLELIC+=($5!~",")} \
# {DBSNP_COUNT+=($3~"rs")} \
# {DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
# END {if (INDEL_BIALLELIC==""&&INDEL_COUNT=="") print "'$SM_TAG'","0","NaN","0","NaN"; \
# else if (INDEL_BIALLELIC==0&&INDEL_COUNT>=1) print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN"; \
# else print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_INDEL_METRICS.TXT"

# ################################################################################
# ##### BASIC METRICS FOR MIXED VARIANT TYPES ON BAIT ############################
# ################################################################################
# ##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT MIXED VARIANT ####################
# ##### THIS IS THE HEADER #######################################################
# ##### "SM_TAG""\t""COUNT_MIXED_ON_BAIT""\t""PERCENT_MIXED_ON_BAIT_SNP138"} ##### 
# ################################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_BAIT/$SM_TAG".MIXED.ON_BAIT.PASS.vcf.gz" \
# | awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
# END {if (MIXED_COUNT!="") print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
# else print "'$SM_TAG'","0","NaN"}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_MIXED_METRICS.TXT"

# ###################################################################################
# ##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET MIXED VARIANT #####################
# ###################################################################################
# ##### THIS IS THE HEADER ##########################################################
# ##### "SM_TAG""\t""COUNT_MIXED_ON_TARGET""\t""PERCENT_MIXED_ON_TARGET_SNP138" ##### 
# ###################################################################################

# zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_TARGET/$SM_TAG".MIXED.ON_TARGET.PASS.vcf.gz" \
# | awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
# END {if (MIXED_COUNT!="") print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
# else print "'$SM_TAG'","0","NaN"}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_MIXED_METRICS.TXT"