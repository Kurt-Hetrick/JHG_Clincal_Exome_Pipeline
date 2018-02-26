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

SAMTOOLS_DIR=$1
CORE_PATH=$2
DATAMASH_DIR=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6
FATHER=$7
MOTHER=$8
GENDER=$9
AFFECTED=${10}

#########################################################
##### Grabbing the BAM header (for RG ID,PU,LB,etc) #####
#########################################################

##############################################################################################################################################
############################################################## OLD ###########################################################################
##############################################################################################################################################
# $SAMTOOLS_DIR/samtools view -H \ ###########################################################################################################
# $CORE_PATH/$PROJECT/BAM/$SM_TAG".bam" \ ####################################################################################################
# | grep ^@RG \ ##############################################################################################################################
# | awk '{split($5,SMtag,":"); split($2,ID,":"); split($6,PU,":"); split($3,Library,":"); print SMtag[2]"\t"ID[2]"\t"PU[2]"\t"Library[2]}' \ #
# | $DATAMASH_DIR/datamash -s -g 1 collapse 2 collapse 3 unique 4 | sed 's/,/;/g' \ ##########################################################
# | awk 'BEGIN {OFS="\t"} {split($4,Library,"_"); print $0,Library[1],"'$FAMILY'","'$FATHER'","'$MOTHER'","'$GENDER'","'$AFFECTED'"}' \ ######
# | awk 'BEGIN {OFS="\t"} $9=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,"MALE",$10} $9=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"FEMALE",$10}' \ #########
# | awk 'BEGIN {OFS="\t"} $10=="-9" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} $10=="0" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"MISSING"} \ #####
# $10=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"UNAFFECTED"} $10=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"AFFECTED"}' \ #########################
# >> $CORE_PATH/$PROJECT/TEMP/SAMPLE_META.txt ################################################################################################
##############################################################################################################################################

#######################################################################################################
##### THIS IS THE HEADER ##############################################################################
##### "PROJECT","SM_TAG","RG_PU","Library_Name","FAMILY","FATHER","MOTHER","LIMS_SEX","PHENOTYPE" #####
#######################################################################################################

$SAMTOOLS_DIR/samtools view -H \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/BAM/$SM_TAG".bam" \
| grep ^@RG \
| awk '{split($5,SMtag,":"); split($6,PU,":"); split($3,Library,":"); print SMtag[2]"\t"PU[2]"\t"Library[2]}' \
| $DATAMASH_DIR/datamash -s -g 1 collapse 2 unique 3 | sed 's/,/;/g' \
| awk 'BEGIN {OFS="\t"} {print $0,"'$PROJECT'","'$FAMILY'","'$FATHER'","'$MOTHER'","'$GENDER'","'$AFFECTED'"}' \
| awk 'BEGIN {OFS="\t"} $8=="1" {print $4,$1,$2,$3,$5,$6,$7,"MALE",$9} $8=="2" {print $4,$1,$2,$3,$5,$6,$7,"FEMALE",$9} $8!="1"&&$8!="2" {print $4,$1,$2,$3,$5,$6,$7,"UNKNOWN",$9}' \
| awk 'BEGIN {OFS="\t"} $9=="-9" {print $1,$2,$3,$4,$5,$6,$7,$8,"MISSING"} $9=="0" {print $1,$2,$3,$4,$5,$6,$7,$8,"MISSING"} \
$9=="1" {print $1,$2,$3,$4,$5,$6,$7,$8,"UNAFFECTED"} $9=="2" {print $1,$2,$3,$4,$5,$6,$7,$8,"AFFECTED"}' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_SAMPLE_META.txt"

########################################################
##### GENDER CHECK FROM ANEUPLOIDY CHECK ###############
########################################################
##### THIS IS THE HEADER ###############################
##### SM_TAG,X_AVG_DP,X_NORM_DP,Y_AVG_DP,Y_NORM_DP #####
########################################################

awk 'BEGIN {OFS="\t"} $2=="X"&&$3=="whole" {print $6,$7} $2=="Y"&&$3=="whole" {print $6,$7}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK/$SM_TAG".chrom_count_report.txt" \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_GENDER_CHECK.TXT"

# # GRABBING CONCORDANCE. MULTI-SAMPLE # This will be just for the validation runs.
# 
# ls $CORE_PATH/$PROJECT/REPORTS/CONCORDANCE_MS/*_concordance.csv \
# | awk '{print "awk","1",$0}' \
# | bash \
# | sort -r \
# | uniq \
# | awk 'NR>1' \
# | sort \
# | sed 's/,/\t/g' \
# | awk 'BEGIN {print "SM_TAG","COUNT_DISC_HOM","COUNT_CONC_HOM","PERCENT_CONC_HOM",\
# "COUNT_DISC_HET","COUNT_CONC_HET","PERCENT_CONC_HET",\
# "PERCENT_TOTAL_CONC","COUNT_HET_BEADCHIP","SENSITIVITY_2_HET"} \
# {print $1,$5,$6,$7,$2,$3,$4,$8,$9,$10}' \
# | sed 's/ /\t/g' \
# >| $CORE_PATH/$PROJECT/TEMP/CONCORDANCE_MS.txt
# 

##############################################################################################################################################
##### VERIFY BAM ID ##########################################################################################################################
##############################################################################################################################################
##### THIS IS THE HEADER #####################################################################################################################
##### "SM_TAG","VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0","VERIFYBAM_DIFF_LK0_LK1","VERIFYBAM_AVG_DP" #####
##############################################################################################################################################

awk 'BEGIN {OFS="\t"} NR>1 {print $1,$7,$4,$8,$9,($9-$8),$6}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/VERIFYBAMID/$SM_TAG".selfSM" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_VERFIY_BAM_ID.TXT"

#############################################################################################
##### INSERT SIZE ###########################################################################
#############################################################################################
##### THIS IS THE HEADER ####################################################################
##### "SM_TAG","MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE" #####
#############################################################################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$1,$5,$6}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_INSERT_SIZE_METRICS.TXT"

##########################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 1 ###############################
##########################################################################
##### THIS THE HEADER ####################################################
##### "SM_TAG","PCT_PF_READS_ALIGNED_R1","PF_HQ_ALIGNED_READS_R1" ########
##### "PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1" #####
##### "PCT_READS_ALIGNED_IN_PAIRS_R1","PCT_ADAPTER_R1" ###################
##########################################################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$7,$9,$13,$14,$15,$18,$22}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_ALIGNMENT_SUMMARY_READ_1_METRICS.TXT"

##########################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 2 ###############################
##########################################################################
##### THIS THE HEADER ####################################################
##### "SM_TAG","PCT_PF_READS_ALIGNED_R2","PF_HQ_ALIGNED_READS_R2" ########
##### "PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2" #####
##### "PCT_READS_ALIGNED_IN_PAIRS_R2","PCT_ADAPTER_R2" ###################
##########################################################################

awk 'BEGIN {OFS="\t"} NR==9 {print "'$SM_TAG'",$7,$9,$13,$14,$15,$18,$22}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_ALIGNMENT_SUMMARY_READ_2_METRICS.TXT"

#######################################################################################
##### ALIGNMENT SUMMARY METRICS FOR PAIR ##############################################
#######################################################################################
##### THIS THE HEADER #################################################################
##### "SM_TAG","TOTAL_READS","RAW_GIGS","PCT_PF_READS_ALIGNED_PAIR" ###################
##### "PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR" ############
##### "PCT_READS_ALIGNED_IN_PAIRS_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR" #####
#######################################################################################

awk 'BEGIN {OFS="\t"} NR==10 {print "'$SM_TAG'",$2,($2*$16/1000000000),$7,$13,$14,$15,$18,$20,$21}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_ALIGNMENT_SUMMARY_READ_PAIR_METRICS.TXT"

###################################################################################################################
##### MARK DUPLICATES REPORT ######################################################################################
###################################################################################################################
##### THIS IS THE HEADER ##########################################################################################
##### "SM_TAG","UNMAPPED_READS","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE" #####
###################################################################################################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$4,$7,$8,$9}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_MARK_DUPLICATES_METRICS.TXT"

##########################################################################################################################
##### HYBRIDIZATION SELECTION REPORT #####################################################################################
##########################################################################################################################
##### THIS IS THE HEADER #################################################################################################
##### "SM_TAG","GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY","PCT_PF_UQ_READS_ALIGNED" ###############################
##### "PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE" #####
##### "ZERO_CVG_TARGETS_PCT","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET" #######################
##### "PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X" ########################
##### "AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q"} #######################################################
##########################################################################################################################

awk 'BEGIN {OFS="\t"} NR==8 {print "'$SM_TAG'",$2,$3,$4,$12,($14/1000000000),$19,$22,$23,$24,$28,$30,$31,$32,$33,\
$38,$39,$40,$41,$50,$51,$52,$53}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_HYB_SELECTION.TXT"

##############################################
##### BAIT BIAS REPORT FOR Cref and Gref #####
##############################################
##### THIS IS THE HEADER #####################
##### SM_TAG,Cref_Q,Gref_Q ###################
##############################################

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt" \
| sed '/^$/d' \
| awk 'BEGIN {OFS="\t"} $12=="Cref"||$12=="Gref"  {print $5}' \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_BIAS.TXT"

############################################################
##### PRE-ADAPTER BIAS REPORT FOR Deamination and OxoG #####
############################################################
##### THIS IS THE HEADER ###################################
##### SM_TAG,Deamination_Q,OxoG_Q ##########################
############################################################

grep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" \
| sed '/^$/d' \
| awk 'BEGIN {OFS="\t"} $12=="Deamination"||$12=="OxoG"  {print $5}' \
| paste - - \
| awk 'BEGIN {OFS="\t"} {print "'$SM_TAG'",$0}' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_PRE_ADAPTER.TXT"

###########################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT SNVS ########################
###########################################################################
##### THIS IS THE HEADER ##################################################
##### "SM_TAG""\t""COUNT_SNV_ON_BAIT""\t""PERCENT_SNV_ON_BAIT_SNP138" ##### 
###########################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_BAIT/$SM_TAG".SNV.ON_BAIT.PASS.vcf.gz" \
| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (SNV_COUNT!="") {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_SNV_METRICS.TXT"

###############################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET SNVS ##########################
###############################################################################
##### THIS IS THE HEADER ######################################################
##### "SM_TAG""\t""COUNT_SNV_ON_TARGET""\t""PERCENT_SNV_ON_TARGET_SNP138" ##### 
###############################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/FILTERED_ON_TARGET/$SM_TAG".SNV.ON_TARGET.PASS.vcf.gz" \
| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (SNV_COUNT!="") {print "'$SM_TAG'",SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
else {print "'$SM_TAG'","0","NaN"}}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_SNV_METRICS.TXT"

##############################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, ALL ###############
##############################################################
##### THIS IS THE HEADER #####################################
##### "SM_TAG""\t""ALL_TI_TV_COUNT""\t""ALL_TI_TV_RATIO" #####
##############################################################

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_All_.titv.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_ALL.TXT"

##################################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, KNOWN #################
##################################################################
##### THIS IS THE HEADER #########################################
##### "SM_TAG""\t""KNOWN_TI_TV_COUNT""\t""KNOWN_TI_TV_RATIO" #####
##################################################################

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Known_.titv.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_KNOWN.TXT"

##################################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, NOVEL #################
##################################################################
##### THIS IS THE HEADER #########################################
##### "SM_TAG""\t""NOVEL_TI_TV_COUNT""\t""NOVEL_TI_TV_RATIO" #####
##################################################################

awk 'BEGIN {OFS="\t"} END {if ($2!="") {print "'$SM_TAG'",$2,$6} \
else {print "'$SM_TAG'","0","NaN"}}' \
$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/TI_TV/$SM_TAG"_Novel_.titv.txt" \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TITV_NOVEL.TXT"

######################################################################################################################################
##### INDEL METRICS ON BAIT ##########################################################################################################
######################################################################################################################################
##### THIS IS THE HEADER #############################################################################################################
##### "SM_TAG","COUNT_ALL_INDEL_BAIT","ALL_INDEL_BAIT_PCT_SNP138","COUNT_BIALLELIC_INDEL_BAIT","BIALLELIC_INDEL_BAIT_PCT_SNP138" #####
######################################################################################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_BAIT/$SM_TAG".INDEL.ON_BAIT.PASS.vcf.gz" \
| awk '{INDEL_COUNT++NR} \
{INDEL_BIALLELIC+=($5!~",")} \
{DBSNP_COUNT+=($3~"rs")} \
{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
END {if (INDEL_BIALLELIC==""&&INDEL_COUNT=="") print "'$SM_TAG'","0","NaN","0","NaN"; \
else if (INDEL_BIALLELIC==0&&INDEL_COUNT>=1) print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN"; \
else print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_INDEL_METRICS.TXT"

##############################################################################################################################################
##### INDEL METRICS ON TARGET ################################################################################################################
##############################################################################################################################################
##### THIS IS THE HEADER #####################################################################################################################
##### "SM_TAG","COUNT_ALL_INDEL_TARGET","ALL_INDEL_TARGET_PCT_SNP138","COUNT_BIALLELIC_INDEL_TARGET","BIALLELIC_INDEL_TARGET_PCT_SNP138" #####
##############################################################################################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/FILTERED_ON_TARGET/$SM_TAG".INDEL.ON_TARGET.PASS.vcf.gz" \
| awk '{INDEL_COUNT++NR} \
{INDEL_BIALLELIC+=($5!~",")} \
{DBSNP_COUNT+=($3~"rs")} \
{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
END {if (INDEL_BIALLELIC==""&&INDEL_COUNT=="") print "'$SM_TAG'","0","NaN","0","NaN"; \
else if (INDEL_BIALLELIC==0&&INDEL_COUNT>=1) print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN"; \
else print "'$SM_TAG'",INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_INDEL_METRICS.TXT"

################################################################################
##### BASIC METRICS FOR MIXED VARIANT TYPES ON BAIT ############################
################################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT MIXED VARIANT ####################
##### THIS IS THE HEADER #######################################################
##### "SM_TAG""\t""COUNT_MIXED_ON_BAIT""\t""PERCENT_MIXED_ON_BAIT_SNP138"} ##### 
################################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_BAIT/$SM_TAG".MIXED.ON_BAIT.PASS.vcf.gz" \
| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (MIXED_COUNT!="") print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
else print "'$SM_TAG'","0","NaN"}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_BAIT_MIXED_METRICS.TXT"

###################################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET MIXED VARIANT #####################
###################################################################################
##### THIS IS THE HEADER ##########################################################
##### "SM_TAG""\t""COUNT_MIXED_ON_TARGET""\t""PERCENT_MIXED_ON_TARGET_SNP138" ##### 
###################################################################################

zgrep -v "^#" $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/FILTERED_ON_TARGET/$SM_TAG".MIXED.ON_TARGET.PASS.vcf.gz" \
| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
END {if (MIXED_COUNT!="") print "'$SM_TAG'",MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
else print "'$SM_TAG'","0","NaN"}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$FAMILY"_TARGET_MIXED_METRICS.TXT"
