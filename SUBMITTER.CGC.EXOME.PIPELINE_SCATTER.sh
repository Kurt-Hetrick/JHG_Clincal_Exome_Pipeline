#!/bin/bash

SAMPLE_SHEET=$1
PED_FILE=$2

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

SCRIPT_DIR="/isilon/cgc/PIPELINES/JHGenomics_CGC_Clinical_Exome/d04bd2a2dc3c00f0ba7fc697750a1eef66985284/scripts"
# The above hash value is the corresponding commit at https://github.com/Kurt-Hetrick/JHGenomics_CGC_Clinical_Exome

CORE_PATH="/isilon/cgc/SS_CRE"
CONTROL_REPO="/isilon/cgc/PIPELINE_FILES"

# PIPELINE PROGRAMS
JAVA_1_6="/isilon/cgc/PROGRAMS/jre1.6.0_25/bin"
JAVA_1_8="/isilon/cgc/PROGRAMS/jdk1.8.0_73/bin"
BWA_DIR="/isilon/cgc/PROGRAMS/bwa-0.7.8"
PICARD_DIR="/isilon/cgc/PROGRAMS/picard-tools-2.1.1"
GATK_DIR="/isilon/cgc/PROGRAMS/GenomeAnalysisTK-3.7"
VERIFY_DIR="/isilon/cgc/PROGRAMS/verifyBamID_20120620/bin/"
TABIX_DIR="/isilon/cgc/PROGRAMS/tabix-0.2.6"
SAMTOOLS_DIR="/isilon/cgc/PROGRAMS/samtools-0.1.18"
DATAMASH_DIR="/isilon/cgc/PROGRAMS/datamash-1.0.6"
BEDTOOLS_DIR="/isilon/cgc/PROGRAMS/bedtools-2.22.0/bin"
VCFTOOLS_DIR="/isilon/cgc/PROGRAMS/vcftools_0.1.12b/bin"
PLINK2_DIR="/isilon/cgc/PROGRAMS/PLINK2"
KING_DIR="/isilon/cgc/PROGRAMS/KING/Linux-king19"
CIDRSEQSUITE_DIR="/isilon/cgc/PROGRAMS/CIDRSeqSuiteSoftware_Version_4_0/"
ANNOVAR_DIR="/isilon/cgc/PROGRAMS/ANNOVAR/2013_09_11"

# PIPELINE FILES
GENE_LIST="/isilon/cgc/PIPELINE_FILES/RefSeqGene.GRCh37.Ready.txt"
VERIFY_VCF="/isilon/cgc/PIPELINE_FILES/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
CODING_BED="/isilon/cgc/PIPELINE_FILES/RefSeq.Unique.GRCh37.FINAL.bed"
CYTOBAND_BED="/isilon/cgc/PIPELINE_FILES/GRCh37.Cytobands.bed"
HAPMAP="/isilon/cgc/PIPELINE_FILES/hapmap_3.3.b37.vcf"
OMNI_1KG="/isilon/cgc/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
HI_CONF_1KG_PHASE1_SNP="/isilon/cgc/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
MILLS_1KG_GOLD_INDEL="/isilon/cgc/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
PHASE3_1KG_AUTOSOMES="/isilon/cgc/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
DBSNP_129="/isilon/cgc/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"

# Can't use /usr/local/sandbox b/c it keeps on getting screwed up since the move to bright

# GENE_LIST="/usr/local/sandbox/cgc/PIPELINE_FILES/RefSeqGene.GRCh37.Ready.txt"
# VERIFY_VCF="/usr/local/sandbox/cgc/PIPELINE_FILES/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
# CODING_BED="/isilon/cgc/PIPELINE_FILES/RefSeq.Unique.GRCh37.FINAL.bed"
# CYTOBAND_BED="/usr/local/sandbox/cgc/PIPELINE_FILES/GRCh37.Cytobands.bed"
# HAPMAP="/usr/local/sandbox/cgc/PIPELINE_FILES/hapmap_3.3.b37.vcf"
# OMNI_1KG="/usr/local/sandbox/cgc/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
# HI_CONF_1KG_PHASE1_SNP="/usr/local/sandbox/cgc/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
# MILLS_1KG_GOLD_INDEL="/usr/local/sandbox/cgc/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
# PHASE3_1KG_AUTOSOMES="/usr/local/sandbox/cgc/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
# DBSNP_129="/usr/local/sandbox/cgc/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"

CONTROL_PED_FILE="/isilon/cgc/PIPELINE_FILES/CGC_CONTROL_SET_UNIQUE.ped"

#################################
##### MAKE A DIRECTORY TREE #####
#################################

##### CREATING A DIRECTORY IN USER'S HOME TO MERGE THE SAMPLE MANIFEST WITH THE PED FILE

mkdir -p ~/CGC_PIPELINE_TEMP

MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
PED_PREFIX=`basename $PED_FILE .ped`

##########################################################

SETUP_PROJECT ()
{
FORMAT_MANIFEST
MERGE_PED_MANIFEST
CREATE_SAMPLE_INFO_ARRAY
MAKE_PROJ_DIR_TREE
# echo "echo Making padded annotated RefSeq coding bed file for $SAMPLE"
# PAD_REFSEQ
# echo "echo Making padded target bed file for $SAMPLE"
# PAD_TARGET
# echo "echo Making everything merged together bait file for $SAMPLE"
# MAKE_BAIT
# echo Project started at `date` >| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/PROJECT_START_END_TIMESTAMP.txt
}

FORMAT_MANIFEST ()
{
sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| sed 's/,/\t/g' \
| sort -k 8,8 \
>| ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt
}

MERGE_PED_MANIFEST ()
{
awk 1 $PED_FILE \
| sed 's/\r//g' \
| sort -k 2,2 \
| join -1 8 -2 2 -e '-'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
}

# MAKE AN ARRAY FOR EACH SAMPLE
	## SAMPLE_INFO_ARRAY[0] = PROJECT
	## SAMPLE_INFO_ARRAY[1] = FAMILY
	## SAMPLE_INFO_ARRAY[2] = SM_TAG
		## SAMPLE = SM_TAG
	## SAMPLE_INFO_ARRAY[3] = BAIT BED FILE
	## SAMPLE_INFO_ARRAY[4] = TARGET_BED_FILE

CREATE_SAMPLE_INFO_ARRAY ()
{
SAMPLE_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$20,$8,$15,$16}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

# PROJECT DIRECTORY TREE CREATOR

MAKE_PROJ_DIR_TREE ()
{
mkdir -p $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/LOGS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/HC_BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/INDEL/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/SNV/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/MIXED/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/GVCF \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_CHR} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/CONCORDANCE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/DEPTH_OF_COVERAGE/{TARGET,REFSEQ_CODING_PLUS_10bp} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/INSERT_SIZE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/${SAMPLE_INFO_ARRAY[2]}/REPORTS/ANEUPLOIDY_CHECK \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/${SAMPLE_INFO_ARRAY[1]}/{LOGS,VCF,RELATEDNESS,PCA} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}_ANNOVAR \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/{VCF_PREP,PLINK,KING} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/{TEMP,FASTQ,REPORTS,LOGS}
}

# PAD THE REFSEQ canonical transcript bed file by 10 bases.
# can make this as an input variable with a default value 10 if i have to ever give more than 0 effs.

PAD_REFSEQ ()
{
awk 1 $CODING_BED \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_CODING.bed"
}

# PAD THE TARGET BED FILE BY 10 BP

PAD_TARGET ()
{
awk 1 ${SAMPLE_INFO_ARRAY[4]} \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_TARGET.bed"
}

# MERGE THE PADDED THE TARGET BED WITH THE BAIT BED FILE

MAKE_BAIT ()
{
cat $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_TARGET.bed" \
${SAMPLE_INFO_ARRAY[3]} \
| sort -k 1,1 -k 2,2n -k 3,3n \
| $BEDTOOLS_DIR/bedtools merge -i - \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}_BAIT.bed
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
SETUP_PROJECT
done

############################################################

# to create the qsub cmd line to submit bwa alignments to the cluster
# handle blank lines
# handle something else too

awk 'BEGIN {FS="\t"} {split($19,INDEL,";");split($8,smtag,"[@]"); \
print "qsub","-N","A.01_BWA_"smtag[1]"_"smtag[2]"_"$2"_"$3"_"$4,\
"-o","'$CORE_PATH'/"$1"/"$20"/"$8"/LOGS/"$8"_"$2"_"$3"_"$4".BWA.log",\
"'$SCRIPT_DIR'""/A.01_BWA.sh",\
"'$BWA_DIR'","'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$20,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12"\n""sleep 3s"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt

# create a hold job id qsub command line based on the number of
# submit merging the bam files created by bwa mem above
# only launch when every lane for a sample is done being processed by bwa mem

awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); print $1,$20,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2,3 collapse 4 collapse 5 \
| awk 'BEGIN {FS="\t"} \
gsub(/,/,",A.01_BWA_"$3"_",$4) \
gsub(/,/,",INPUT=/isilon/cgc/SS_CRE/"$1"/TEMP/",$5) \
{split($3,smtag,"[@]"); print "qsub","-N","B.01_MERGE_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".MERGE.BAM.FILES.log",\
"-hold_jid","A.01_BWA_"smtag[1]"_"smtag[2]"_"$4, \
"'$SCRIPT_DIR'""/B.01_MERGE_SORT_AGGRO.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3,"INPUT=/isilon/cgc/SS_CRE/"$1"/TEMP/"$5"\n""sleep 1s"}'

# Mark duplicates on the bam file above. Create a Mark Duplicates report which goes into the QC report

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","C.01_MARK_DUPLICATES_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","B.01_MERGE_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".MARK_DUPLICATES.log",\
"'$SCRIPT_DIR'""/C.01_MARK_DUPLICATES.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# Generate a list of places that could be potentially realigned.

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
print "qsub","-N","D.01_REALIGNER_TARGET_CREATOR_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","C.01_MARK_DUPLICATES_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".REALIGNER_TARGET_CREATOR.log",\
"'$SCRIPT_DIR'""/D.01_REALIGNER_TARGET_CREATOR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2]"\n""sleep 1s"}'

# With the list generated above walk through the BAM file and realign where necessary
# Write out a new bam file

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
print "qsub","-N","E.01_INDEL_REALIGNER_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","D.01_REALIGNER_TARGET_CREATOR_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".INDEL_REALIGNER.log",\
"'$SCRIPT_DIR'""/E.01_INDEL_REALIGNER.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2]"\n""sleep 1s"}'

# Run Base Quality Score Recalibration

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19,$18,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
print "qsub","-N","F.01_PERFORM_BQSR_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","E.01_INDEL_REALIGNER_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PERFORM_BQSR.log",\
"'$SCRIPT_DIR'""/F.01_PERFORM_BQSR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6,$7"\n""sleep 1s"}'

# write Final Bam file

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","F.01_PERFORM_BQSR_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".FINAL_BAM.log",\
"'$SCRIPT_DIR'""/G.01_FINAL_BAM.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# SCATTER THE HAPLOTYPE CALLER GVCF CREATION USING THE WHERE THE BED INTERSECTS WITH {{1.22},{X,Y}}

CREATE_SAMPLE_INFO_ARRAY_HC ()
{
SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_HAPLOTYPE_CALLER ()
{
echo \
qsub \
-N H.01_HAPLOTYPE_CALLER_${SAMPLE_INFO_ARRAY_HC[0]}_${SAMPLE_INFO_ARRAY_HC[5]}_chr$CHROMOSOME \
-hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_HC[5]}_${SAMPLE_INFO_ARRAY_HC[0]} \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_chr$CHROMOSOME.log \
$SCRIPT_DIR/H.01_HAPLOTYPE_CALLER_SCATTER.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]} \
$CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_HC
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_HAPLOTYPE_CALLER
		echo sleep 1s
		done
	done

################################################################

# GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES INTO A SINGLE SAMPLE GVCF

# GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES INTO A SINGLE SAMPLE GVCF

# BUILD_HOLD_ID_PATH(){
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"H.01_HAPLOTYPE_CALLER_"$PROJECT"_"${SAMPLE_INFO_ARRAY_HC[5]}"_chr"$CHROMOSOME","
#  	done
#  done
# }

BUILD_HOLD_ID_PATH(){
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"H.01_HAPLOTYPE_CALLER_"$PROJECT"_"$SAMPLE"_chr"$CHROMOSOME","
 	done
 done
}

# CREATE_SAMPLE_INFO_ARRAY_HC ()
# {
# SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

CREATE_SAMPLE_INFO_ARRAY_HC ()
{
SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_HAPLOTYPE_CALLER_GATHER ()
{
echo \
qsub \
-N H.01-A.01_HAPLOTYPE_CALLER_GATHER_${SAMPLE_INFO_ARRAY_HC[0]}_$SAMPLE \
${HOLD_ID_PATH} \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_GATHER.log \
$SCRIPT_DIR/H.01-A.01_HAPLOTYPE_CALLER_GATHER.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]}
}

# CALL_HAPLOTYPE_CALLER_GATHER ()
# {
# echo \
# qsub \
# -N H.01-A.01_HAPLOTYPE_CALLER_GATHER_${SAMPLE_INFO_ARRAY_HC[0]}_${SAMPLE_INFO_ARRAY_HC[5]} \
# ${HOLD_ID_PATH} \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_GATHER.log \
# $SCRIPT_DIR/H.01-A.01_HAPLOTYPE_CALLER_GATHER.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]}
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8} $SAMPLE_SHEET | sort | uniq );
#  do
# 	BUILD_HOLD_ID_PATH
# 	CREATE_SAMPLE_INFO_ARRAY_HC
# 	CALL_HAPLOTYPE_CALLER_GATHER
# 	echo sleep 1s
#  done

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]}' $SAMPLE_SHEET | sort | uniq );

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
 do
	BUILD_HOLD_ID_PATH
	CREATE_SAMPLE_INFO_ARRAY_HC
	CALL_HAPLOTYPE_CALLER_GATHER
	echo sleep 1s
 done

# Run POST BQSR TABLE

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
print "qsub","-N","H.02_POST_BQSR_TABLE_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".POST_BQSR_TABLE.log",\
"'$SCRIPT_DIR'""/H.02_POST_BQSR_TABLE.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6"\n""sleep 1s"}'

# Run ANALYZE COVARIATES

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
print "qsub","-N","H.02-A.01_ANALYZE_COVARIATES_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.02_POST_BQSR_TABLE_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANALYZE_COVARIATES.log",\
"'$SCRIPT_DIR'""/H.02-A.01_ANALYZE_COVARIATES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN DOC CODING PLUS 10 BP FLANKS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_CODING_10bpFLANKS.log",\
"'$SCRIPT_DIR'""/H.03_DOC_CODING_10bpFLANKS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN ANEUPLOIDY_CHECK AFTER CODING PLUS 10 BP FLANKS FINISHES

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03-A.01_DOC_CHROM_DEPTH_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANEUPLOIDY_CHECK.log",\
"'$SCRIPT_DIR'""/H.03-A.01_CHROM_DEPTH.sh",\
"'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2,$3"\n""sleep 1s"}'

# RUN FORMATTING PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE.log",\
"'$SCRIPT_DIR'""/H.03-A.02_PER_BASE.sh",\
"'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2,$3"\n""sleep 1s"}'

# RUN FILTERING PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_FILTER.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.01_PER_BASE_FILTERED.sh",\
"'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_BGZIP.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.02_PER_BASE_BGZIP.sh",\
"'$CORE_PATH'","'$TABIX_DIR'",$1,$2,$3"\n""sleep 1s"}'

# TABIX PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.03-A.02-A.02-A.01_PER_BASE_TABIX_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_TABIX.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.02-A.01_PER_BASE_TABIX.sh",\
"'$CORE_PATH'","'$TABIX_DIR'",$1,$2,$3"\n""sleep 1s"}'

# RUN DOC TARGET BED

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.05_DOC_TARGET_BED_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_TARGET_BED.log",\
"'$SCRIPT_DIR'""/H.05_DOC_TARGET_BED.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN COLLECT MULTIPLE METRICS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$18,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.06_COLLECT_MULTIPLE_METRICS_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_MULTIPLE_METRICS.log",\
"'$SCRIPT_DIR'""/H.06_COLLECT_MULTIPLE_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5,$6"\n""sleep 1s"}'

# RUN COLLECT HS METRICS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$16,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.07_COLLECT_HS_METRICS_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_HS_METRICS.log",\
"'$SCRIPT_DIR'""/H.07_COLLECT_HS_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN SELECT VERIFYBAM ID VCF

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.08_SELECT_VERIFYBAMID_VCF_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".SELECT_VERIFYBAMID_VCF.log",\
"'$SCRIPT_DIR'""/H.08_SELECT_VERIFYBAMID_VCF.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VERIFY_VCF'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# RUN VERIFYBAMID

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); \
print "qsub","-N","H.08-A.01_VERIFYBAMID_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.08_SELECT_VERIFYBAMID_VCF_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".VERIFYBAMID.log",\
"'$SCRIPT_DIR'""/H.08-A.01_VERIFYBAMID.sh",\
"'$CORE_PATH'","'$VERIFY_DIR'",$1,$2,$3"\n""sleep 1s"}'

###################################################
### RUN VERIFYBAM ID PER CHROMOSOME - VITO ########
###################################################

CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM ()
{
SAMPLE_INFO_ARRAY_VERIFY_BAM=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$15,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_SELECT_VERIFY_BAM ()
{
echo \
qsub \
-N H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.SELECT_VERIFYBAMID_chr$CHROMOSOME.log \
$SCRIPT_DIR/H.09_SELECT_VERIFYBAMID_VCF_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH $VERIFY_VCF \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[3]} \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[4]} $CHROMOSOME
}

CALL_VERIFYBAMID ()
{
echo \
qsub \
-N H.09-A.01_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-hold_jid H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.VERIFYBAMID_chr$CHROMOSOME.log \
$SCRIPT_DIR/H.09-A.01_VERIFYBAMID_CHR.sh \
$CORE_PATH $VERIFY_DIR \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} \
$CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
	for CHROMOSOME in {1..22}
		do
		CALL_SELECT_VERIFY_BAM
		echo sleep 1s
		CALL_VERIFYBAMID
		echo sleep 1s
	done
done

#####################################################
### JOIN THE PER CHROMOSOME VERIFYBAMID REPORTS #####
#####################################################

BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"H.09-A.01_VERIFYBAMID_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}"_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}"_"chr$CHROMOSOME","
 	done
 done
}

 CAT_VERIFYBAMID_CHR ()
 {
echo \
qsub \
-N H.09-A.01-A.01_JOIN_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
$HOLD_ID_PATH \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.CAT_VERIFYBAMID_CHR.log \
$SCRIPT_DIR/H.09-A.01-A.01_CAT_VERIFYBAMID_CHR.sh \
$CORE_PATH \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}
 }

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
 do
 	CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
	BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR
	CAT_VERIFYBAMID_CHR
	echo sleep 1s
 done

#### JOINT CALLING AND VQSR #### ###VITO###

### CREATE A GVCF ".list" file for each sample

CREATE_FAMILY_INFO_ARRAY ()
{
FAMILY_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {split($8,smtag,"[@]"); print $1,smtag[1]"_"smtag[2],$20,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CREATE_GVCF_LIST ()
{
awk 'BEGIN {OFS="/"} $20=="'$FAMILY'" {print "'$CORE_PATH'",$1,$20,$8,"GVCF",$8".g.vcf.gz"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort \
| uniq \
>| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".gvcf.list"
}

CREATE_FAMILY_SAMPLE_LIST ()
{
awk '$20=="'$FAMILY'" {print $8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort \
| uniq \
>| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".sample.list"
}

BUILD_HOLD_ID_PATH_GENOTYPE_GVCF ()
{
##NEED FULL LIST OF SAMPLES IN FAMILY FOR THE INNER FOR LOOP.  THE ${FAMILY_INFO_ARRAY[1]} WON'T WORK CONSIDERING THAT IS SPECIFIC FOR ONE SAMPLE

for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
do
HOLD_ID_PATH="-hold_jid "
for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
	do
		HOLD_ID_PATH=$HOLD_ID_PATH"H.01-A.01_HAPLOTYPE_CALLER_GATHER_"$PROJECT"_"${FAMILY_INFO_ARRAY[1]}","
	done
done
}

CALL_GENOTYPE_GVCF ()
{
for CHROM in {{1..22},{X,Y}};
do
echo \
qsub \
-N I.01_GENOTYPE_GVCF_SCATTER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}_chr$CHROM \
$HOLD_ID_PATH \
-o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.GENOTYPE_GVCF_chr$CHROM.log \
$SCRIPT_DIR/I.01_GENOTYPE_GVCF_SCATTER.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]} ${FAMILY_INFO_ARRAY[4]} $CHROM $CONTROL_REPO
echo sleep 1s
done
}

for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
do
	 CREATE_FAMILY_INFO_ARRAY
	 CREATE_GVCF_LIST
 CREATE_FAMILY_SAMPLE_LIST
BUILD_HOLD_ID_PATH_GENOTYPE_GVCF
CALL_GENOTYPE_GVCF
done

########################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME GVCF FILES INTO A SINGLE FAMILY GVCF #####
########################################################################################

BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHR in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"I.01_GENOTYPE_GVCF_SCATTER_"$FAMILY"_"$PROJECT"_chr"$CHR","
 	done
 done
}


CALL_GENOTYPE_GVCF_GATHER ()
{
echo \
qsub \
-N I.01-A.01_GENOTYPE_GVCF_GATHER_${FAMILY_INFO_ARRAY[0]}_${FAMILY_INFO_ARRAY[2]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.GENOTYPE_GVCF_GATHER.log \
 $SCRIPT_DIR/I.01-A.01_GENOTYPE_GVCF_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
}


for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
 do
 	# echo $FAMILY
	BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER
	CREATE_FAMILY_INFO_ARRAY
	CALL_GENOTYPE_GVCF_GATHER
	echo sleep 1s
 done

##########################################################
################END VITO##################################
##########################################################

#####################################################################################################
##### Run Variant Recalibrator for the SNP model, this is done in parallel with the INDEL model #####
#####################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","J.01_VARIANT_RECALIBRATOR_SNP_"$2"_"$1,\
"-hold_jid","I.01-A.01_GENOTYPE_GVCF_GATHER_"$1"_"$2,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_SNP.log",\
"'$SCRIPT_DIR'""/J.01_VARIANT_RECALIBRATOR_SNP.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,"'$HAPMAP'","'$OMNI_1KG'","'$HI_CONF_1KG_PHASE1_SNP'""\n""sleep 1s"}'

### Run Variant Recalibrator for the INDEL model, this is done in parallel with the SNP model

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","J.02_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
"-hold_jid","I.01-A.01_GENOTYPE_GVCF_GATHER_"$1"_"$2,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_INDEL.log",\
"'$SCRIPT_DIR'""/J.02_VARIANT_RECALIBRATOR_INDEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,"'$MILLS_1KG_GOLD_INDEL'""\n""sleep 1s"}'

### Run Apply Recalbration with the SNP model to the VCF file

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","K.01_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
"-hold_jid","J.01_VARIANT_RECALIBRATOR_SNP_"$2"_"$1",""J.02_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_SNP.log",\
"'$SCRIPT_DIR'""/K.01_APPLY_RECALIBRATION_SNP.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

### Run Apply Recalibration with the INDEL model to the VCF file.

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","L.01_APPLY_RECALIBRATION_INDEL_"$2"_"$1,\
"-hold_jid","K.01_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_INDEL.log",\
"'$SCRIPT_DIR'""/L.01_APPLY_RECALIBRATION_INDEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

################################################
##### SCATTER GATHER FOR ADDING ANNOTATION #####
################################################

CREATE_FAMILY_INFO_ARRAY ()
{
FAMILY_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $1,$8,$20,$12,$18}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_VARIANT_ANNOTATOR ()
{
echo \
qsub \
-N P.01_VARIANT_ANNOTATOR_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}_$CHROMOSOME \
-hold_jid L.01_APPLY_RECALIBRATION_INDEL_${FAMILY_INFO_ARRAY[2]}"_"${FAMILY_INFO_ARRAY[0]} \
-o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_ANNOTATOR_$CHROMOSOME.log \
$SCRIPT_DIR/P.01_VARIANT_ANNOTATOR_SCATTER.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH $PED_FILE \
${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]} $CHROMOSOME $PHASE3_1KG_AUTOSOMES
}

for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
do
CREATE_FAMILY_INFO_ARRAY
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_VARIANT_ANNOTATOR
		echo sleep 1s
	done
done

##############################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME ANNOTATED VCF FILES INTO A SINGLE VCF FILE #####
##############################################################################################

BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"P.01_VARIANT_ANNOTATOR_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

CALL_VARIANT_ANNOTATOR_GATHER ()
{
echo \
qsub \
-N P.01-A.01_VARIANT_ANNOTATOR_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.ADD_MORE_ANNOTATION_GATHER.log \
 $SCRIPT_DIR/P.01-A.01_VARIANT_ANNOTATOR_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
}


for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
 do
	BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION
	CREATE_FAMILY_INFO_ARRAY
	CALL_VARIANT_ANNOTATOR_GATHER
	echo sleep 1s
 done

############################################################################################################
##### DO PER CHROMOSOME VARIANT TO TABLE FOR COHORT ########################################################
############################################################################################################

CREATE_FAMILY_ONLY_ARRAY ()
{
FAMILY_ONLY_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $1,$20,$12,$18}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_VARIANT_TO_TABLE_COHORT_ALL_SITES ()
{
echo \
qsub \
-N P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-hold_jid P.01_VARIANT_ANNOTATOR_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.VARIANT_TO_TABLE_COHORT_ALL_SITES_$CHROMOSOME.log \
$SCRIPT_DIR/P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
}

for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
do
CREATE_FAMILY_ONLY_ARRAY
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_VARIANT_TO_TABLE_COHORT_ALL_SITES
		echo sleep 1s
		done
	done

################################################################################################################
##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR COHORT ########################################################
################################################################################################################

BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_COHORT_GATHER ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

CALL_VARIANT_TO_TABLE_COHORT_GATHER ()
{
echo \
qsub \
-N T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER.log \
 $SCRIPT_DIR/T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
}

for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
 do
	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_COHORT_GATHER
	CREATE_FAMILY_INFO_ARRAY
	CALL_VARIANT_TO_TABLE_COHORT_GATHER
	echo sleep 1s
 done

##############################################################################################################
## BGZIP INITIAL JOINT CALLED VCF TABLE ######################################################################
##############################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$2"_"$1,\
"-hold_jid","T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

##############################################################################################################
## TABIX INDEX INITIAL JOINT CALLED VCF TABLE ################################################################
##############################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","T.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES_"$2"_"$1,\
"-hold_jid","T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'


#################################################################################
########### RUNNING FILTER TO FAMILY ALL SITES BY CHROMOSOME ####################
#################################################################################

CALL_FILTER_TO_FAMILY_ALL_SITES ()
{
echo \
qsub \
-N P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-hold_jid P.01_VARIANT_ANNOTATOR_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.FILTER_TO_FAMILY_ALL_SITES_$CHROMOSOME.log \
$SCRIPT_DIR/P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
}

for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
do
CREATE_FAMILY_ONLY_ARRAY
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_FILTER_TO_FAMILY_ALL_SITES
		echo sleep 1s
		done
	done
	
#####################################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME FILTER TO FAMILY VCF FILES INTO A SINGLE VCF FILE #####
#####################################################################################################

BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

CALL_FILTER_TO_FAMILY_VCF_GATHER ()
{
echo \
qsub \
-N T.03-1_FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.FILTER_TO_FAMILY_ALL_SITES_GATHER.log \
 $SCRIPT_DIR/T.03-1_FILTER_TO_FAMILY_ALL_SITES_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
}

for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
 do
	BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF
	CREATE_FAMILY_INFO_ARRAY
	CALL_FILTER_TO_FAMILY_VCF_GATHER
	echo sleep 1s
 done

############################################################################################################
##### DO PER CHROMOSOME VARIANT TO TABLE FOR FAMILY ########################################################
############################################################################################################

CALL_VARIANT_TO_TABLE_FAMILY_ALL_SITES ()
{
echo \
qsub \
-N T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-hold_jid P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
-o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.VARIANT_TO_TABLE_FAMILY_ALL_SITES_$CHROMOSOME.log \
$SCRIPT_DIR/T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
}

for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
do
CREATE_FAMILY_ONLY_ARRAY
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_VARIANT_TO_TABLE_FAMILY_ALL_SITES
		echo sleep 1s
		done
	done

################################################################################################################
##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR FAMILY ########################################################
################################################################################################################

BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_FAMILY_GATHER ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}


CALL_VARIANT_TO_TABLE_FAMILY_GATHER ()
{
echo \
qsub \
-N T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_TO_TABLE_ALL_SITES_GATHER.log \
 $SCRIPT_DIR/T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
}

for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
 do
	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_FAMILY_GATHER
	CREATE_FAMILY_INFO_ARRAY
	CALL_VARIANT_TO_TABLE_FAMILY_GATHER
	echo sleep 1s
 done

##############################################################################################################
## BGZIP FAMILY ONLY VCF TABLE ###############################################################################
##############################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES_"$2"_"$1,\
"-hold_jid","T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

##############################################################################################################
## TABIX INDEX FAMILY ONLY VCF TABLE #########################################################################
##############################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","T.03-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES_"$2"_"$1,\
"-hold_jid","T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.03-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

#################################################################################
########### RUNNING FILTER TO SAMPLE ALL SITES BY CHROMOSOME ####################
#################################################################################

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$8,$20,$12,$18,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_2
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_FILTER_TO_SAMPLE_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done

CREATE_SAMPLE_INFO_ARRAY_2 ()
{
SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_FILTER_TO_SAMPLE_ALL_SITES ()
{
echo \
qsub \
-N P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-hold_jid P.01_VARIANT_ANNOTATOR_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_$CHROMOSOME.log \
$SCRIPT_DIR/P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_2
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_FILTER_TO_SAMPLE_ALL_SITES
		echo sleep 1s
		done
	done

#####################################################################################################
##### GATHER UP THE PER SAMPLE PER CHROMOSOME FILTER TO SAMPLE VCF FILES INTO A SINGLE VCF FILE #####
#####################################################################################################

BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

# CALL_FILTER_TO_SAMPLE_VCF_GATHER ()
# {
# echo \
# qsub \
# -N T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[2]}/${SAMPLE_INFO_ARRAY_2[1]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[3]}
# }

# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}'

CALL_FILTER_TO_SAMPLE_VCF_GATHER ()
{
echo \
qsub \
-N T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_GATHER.log \
 $SCRIPT_DIR/T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
}

# for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
#  	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF
# 	CREATE_SAMPLE_INFO_ARRAY_2
# 	CALL_FILTER_TO_SAMPLE_VCF_GATHER
# 	echo sleep 1s
#  done

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
 do
 	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF
	CREATE_SAMPLE_INFO_ARRAY_2
	CALL_FILTER_TO_SAMPLE_VCF_GATHER
	echo sleep 1s
 done

############################################################################################################
##### DO PER CHROMOSOME VARIANT TO TABLE FOR SAMPLE ########################################################
############################################################################################################

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

CALL_VARIANT_TO_TABLE_SAMPLE_ALL_SITES ()
{
echo \
qsub \
-N T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-hold_jid P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.VARIANT_TO_TABLE_SAMPLE_ALL_SITES_$CHROMOSOME.log \
$SCRIPT_DIR/T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_2
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_VARIANT_TO_TABLE_SAMPLE_ALL_SITES
		echo sleep 1s
		done
	done

################################################################################################################
##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR SAMPLE ########################################################
################################################################################################################

BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_SAMPLE_GATHER ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

CALL_VARIANT_TO_TABLE_SAMPLE_GATHER ()
{
echo \
qsub \
-N T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER.log \
 $SCRIPT_DIR/T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
 do
	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_SAMPLE_GATHER
	CREATE_SAMPLE_INFO_ARRAY_2
	CALL_VARIANT_TO_TABLE_SAMPLE_GATHER
	echo sleep 1s
 done

#################################################################################################################
## ## BGZIP SAMPLE ONLY VCF TABLE ###############################################################################
#################################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

#################################################################################################################
## ## TABIX INDEX SAMPLE ONLY VCF TABLE #########################################################################
#################################################################################################################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","T.06-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.log",\
"'$SCRIPT_DIR'""/T.06-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

###########################################################################################
########### RUNNING FILTER TO SAMPLE ALL SITES BY CHROMOSOME ON TARGET ####################
###########################################################################################

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$8,$20,$12,$18}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

CREATE_SAMPLE_INFO_ARRAY_2 ()
{
SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
print $1,$20,$8,$12,$18,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_FILTER_TO_SAMPLE_ALL_SITES_ON_TARGET ()
{
echo \
qsub \
-N P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-hold_jid P.01_VARIANT_ANNOTATOR_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_TARGET_$CHROMOSOME.log \
$SCRIPT_DIR/P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH \
${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_2
	for CHROMOSOME in {{1..22},{X,Y}}
		do
		CALL_FILTER_TO_SAMPLE_ALL_SITES_ON_TARGET
		echo sleep 1s
		done
	done

###############################################################################################################
##### GATHER UP THE PER SAMPLE PER CHROMOSOME FILTER TO SAMPLE VCF FILES ON TARGET INTO A SINGLE VCF FILE #####
###############################################################################################################

BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF_TARGET ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
 	done
 done
}

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$18,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

CALL_FILTER_TO_SAMPLE_VCF_TARGET_GATHER ()
{
echo \
qsub \
-N T.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]} \
 ${HOLD_ID_PATH} \
 -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER.log \
 $SCRIPT_DIR/T.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER.sh \
 $JAVA_1_8 $GATK_DIR $CORE_PATH \
 ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
 do
 	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF_TARGET
	CREATE_SAMPLE_INFO_ARRAY_2
	CALL_FILTER_TO_SAMPLE_VCF_TARGET_GATHER
	echo sleep 1s
 done

###############################
##### DOING VCF BREAKOUTS #####
###############################

### SUBSETTING FROM COHORT (FAMILY PLUS CONTROL SET) VCF ###

# FILTER TO JUST VARIANT SITES
# I think Molly might like this output, but if not, then don't have to generate it.

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.01_FILTER_COHORT_VARIANT_ONLY_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY.log",\
"'$SCRIPT_DIR'""/S.01_FILTER_COHORT_VARIANT_ONLY.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# FILTER TO JUST PASSING VARIANT SITES
# I think statgen is using this for some of their programs
# If not needed then don't generate

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.02_FILTER_COHORT_VARIANT_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.02_FILTER_COHORT_VARIANT_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# FILTER TO JUST PASSING BIALLELIC SNV SITES
# TEMPORARY FILE USED FOR PCA AND RELATEDNESS

awk 'BEGIN {OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{print "qsub","-N","S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC.log",\
"'$SCRIPT_DIR'""/S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# RUN HUAS WORKFLOW FOR PCA AND RELATEDNESS

awk 'BEGIN {OFS="\t"} {print $1,$20,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{print "qsub","-N","S.03-A.01_PCA_RELATEDNESS_"$2"_"$1,\
"-hold_jid","S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".PCA_RELATEDNESS.log",\
"'$SCRIPT_DIR'""/S.03-A.01_PCA_RELATEDNESS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VCFTOOLS_DIR'","'$PLINK2_DIR'","'$KING_DIR'",$1,$2,$3,"'$PED_FILE'","'$CONTROL_PED_FILE'""\n""sleep 1s"}'

#################################
### SUBSETTING TO SAMPLE VCFS ###
#################################

## SUBSET TO SAMPLE VARIANTS ONLY ON BAIT

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.07_FILTER_TO_SAMPLE_VARIANTS_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".FILTER_TO_VARIANTS.log",\
"'$SCRIPT_DIR'""/S.07_FILTER_TO_SAMPLE_VARIANTS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING SNVS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.09_FILTER_TO_SAMPLE_SNV_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO SAMPLE PASSING INDELS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.10_FILTER_TO_INDEL_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.10_FILTER_TO_SAMPLE_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO SAMPLE PASSING MIXED

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.11_FILTER_TO_MIXED_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.11_FILTER_TO_SAMPLE_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO TARGET SNV ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO TARGET INDEL ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO TARGET MIXED ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

## SUBSET TO SAMPLE VARIANTS ONLY ON TARGET

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".FILTER_TO_VARIANTS_TARGET.log",\
"'$SCRIPT_DIR'""/S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'


####################
### TITV SECTION ###
####################

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF.log",\
"'$SCRIPT_DIR'""/S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND OVERLAP WITH DBSNP 129

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_KNOWN.log",\
"'$SCRIPT_DIR'""/S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,"'$DBSNP_129'""\n""sleep 1s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND DO NOT OVERLAP WITH DBSNP 129

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_NOVEL.log",\
"'$SCRIPT_DIR'""/S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,"'$DBSNP_129'""\n""sleep 1s"}'

### RUN TITV FOR THE PASSING SNVS THAT FALL IN UCSC CODING REGIONS THAT TOUCH EITHER THE BED OR TARGET FILE

## ALL SNVS TITV

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.01-A.01_TITV_ALL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_ALL.log",\
"'$SCRIPT_DIR'""/S.09-A.01-A.01_TITV_ALL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

## ALL KNOWN SNVS TITV

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.02-A.01_TITV_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_KNOWN.log",\
"'$SCRIPT_DIR'""/S.09-A.02-A.01_TITV_KNOWN.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

## ALL NOVEL SNVS TITV

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.03-A.01_TITV_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_NOVEL.log",\
"'$SCRIPT_DIR'""/S.09-A.03-A.01_TITV_NOVEL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

###################
##### ANNOVAR #####
###################

## RUN ANNOVAR

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 -k 3,3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16-A.01_RUN_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET_"smtag[1]"_"smtag[2]"_"$1,\
"-pe slots 5",\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_ANNOVAR.log",\
"'$SCRIPT_DIR'""/S.16-A.01_RUN_ANNOVAR.sh",\
"'$JAVA_1_6'","'$CIDRSEQSUITE_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## REFORMAT ANNOVAR

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 -k 3,3 \
| uniq \
| awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16-A.01-A.01_REFORMAT_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-hold_jid","S.16-A.01_RUN_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".REFORMAT_ANNOVAR.log",\
"'$SCRIPT_DIR'""/S.16-A.01-A.01_REFORMAT_ANNOVAR.sh",\
"'$ANNOVAR_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

######### FINISH UP #################

### QC REPORT PREP ###

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$21,$22,$23,$24}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk 'BEGIN {FS="\t"}
{split($3,smtag,"[@]"); print "qsub","-N","X.01-QC_REPORT_PREP_"$1"_"smtag[1]"_"smtag[2],\
"-hold_jid","S.16-A.01-A.01_REFORMAT_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".QC_REPORT_PREP.log",\
"'$SCRIPT_DIR'""/X.01-QC_REPORT_PREP.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'","'$DATAMASH_DIR'",$1,$2,$3,$4,$5,$6,$7"\n""sleep 1s"}'

### END PROJECT TASKS ###

awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); print $1,smtag[1]"_"smtag[2]}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1 collapse 2 \
| awk 'BEGIN {FS="\t"}
gsub (/,/,",X.01-QC_REPORT_PREP_"$1"_",$2) \
{print "qsub","-N","X.01-X.01-END_PROJECT_TASKS_"$1,\
"-hold_jid","X.01-QC_REPORT_PREP_"$1"_"$2,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".END_PROJECT_TASKS.log",\
"'$SCRIPT_DIR'""/X.01-X.01-END_PROJECT_TASKS.sh",\
"'$CORE_PATH'","'$DATAMASH_DIR'",$1"\n""sleep 1s"}'

### kEY FOR BLAH ###
#
#      1  CGC_160212_HJLWVBCXX_CGCDev1_TEST
#      2  HJLWVBCXX
#      3  1
#      4  ATGCCTAA
#      5  ILLUMINA
#      6  A01_NA12878
#      7  2/12/2016
#      8  NA12878
#      9  CGC
#     10  HiSeq2500_RapidRun
#     11  HJLWVBCXX_1_ATGCCTAA_A01_NA12878
#     12  /isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta
#     13  MBS
#     14  -2
#     15  /isilon/sequencing/data/Work/BED/Production_BED_files/TsTv_BED_File_Agilent_ClinicalExome_S06588914_OnExon_merged_021015_noCHR.bed
#     16  /isilon/sequencing/data/Work/BED/Production_BED_files/ALLBED_BED_File_Agilent_ClinicalExome_S06588914_ALLBed_merged_021015_noCHR.bed
#     17  /isilon/sequencing/data/Work/BED/Production_BED_files/Targets_BED_File_Agilent_ClinicalExome_S06588914_OnTarget_merged_noCHR_013015.bed
#     18  /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.vcf
#     19  /isilon/sequencing/GATK_resource_bundle/2.2/b37/1000G_phase1.indels.b37.vcf;/isilon/sequencing/GATK_resource_bundle/2.2/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
#     20  XC01463
#     21  NA12891
#     22  NA12892
#     23  2
#     24  2

#######


###### SAMPLE MANIFEST KEY...NOT SURE WHAT I AM GOING TO END UP DOING HERE ######

# PROJECT=$1 # the Seq Proj folder name. 1st column in sample manifest
# FLOWCELL=$2 # flowcell that sample read group was performed on. 2nd column of sample manifest
# LANE=$3 # lane of flowcell that sample read group was performed on. 3rd column of the sample manifest
# INDEX=$4 # sample barcode. 4th column of the sample manifest
# PLATFORM=$5 # type of sequencing chemistry matching SAM specification. 5th column of the sample manifest.
# LIBRARY_NAME=$6 # library group of the sample read group.
# 								# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not
# 								# 6th column of the sample manifest
# RUN_DATE=$7 # should be the run set up date to match the seq run folder name, but it has been arbitrarily populated. field X of manifest.
# SM_TAG=$8 # sample ID. sample name for all files, etc. field X of manifest
# CENTER=$9 # the center/funding mechanism. field X of manifest.
# DESCRIPTION=${10} # Generally we use to denote the sequencer setting (e.g. rapid run). field X of manifest.
# REF_GENOME=${11} # the reference genome used in the analysis pipeline. field X of manifest.
# TI_TV_BED=${12} # populated from sample manifest. where ucsc coding exons overlap with bait and target bed files
# BAIT_BED=${13} # populated from sample manifest. a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
# 								# Used for limited where to run base quality score recalibration on where to create gvcf files.
# TARGET_BED=${14} # populated from sample manifest. bed file acquired from manufacturer of their targets. field X of sample manifest.
# DBSNP=${15} # populated from sample manifest. used to annotate ID field in VCF file. masking in base call quality score recalibration.
# KNOWN_INDEL_1=${16} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
# KNOWN_INDEL_2=${17} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
#
# RIS_ID=${SM_TAG%@*} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
# BARCODE_2D=${SM_TAG#*@} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
#
####################################################################################
