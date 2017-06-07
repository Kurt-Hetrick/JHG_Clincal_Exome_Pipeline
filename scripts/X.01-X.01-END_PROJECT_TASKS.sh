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
DATAMASH=$2

PROJECT=$3

# Sorting concatenated sample meta by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*SAMPLE_META.txt \
| sort -k 2,2 \
| uniq \
| awk 'BEGIN {print "PROJECT","SM_TAG","RG_PU","LIBRARY","FAMILY","FATHER","MOTHER","LIMS_SEX","PHENOTYPE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/SAMPLE_META_HEADER.txt

# Sorting concatenated gender check report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*GENDER_CHECK.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","X_AVG_DP","X_NORM_DP","Y_AVG_DP","Y_NORM_DP"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK_HEADER.TXT

# Sorting concatenated verify bam ID report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*VERFIY_BAM_ID.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0","VERIFYBAM_DIFF_LK0_LK1","VERIFYBAM_AVG_DP"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID_HEADER.TXT

# Sorting concatenated insert size report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*INSERT_SIZE_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS_HEADER.TXT

# Sorting concatenated read 1 alignment summary metric report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*ALIGNMENT_SUMMARY_READ_1_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","PCT_PF_READS_ALIGNED_R1","PF_HQ_ALIGNED_READS_R1",\
"PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1",\
"PCT_READS_ALIGNED_IN_PAIRS_R1","PCT_ADAPTER_R1"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS_HEADER.TXT

# Sorting concatenated read 2 alignment summary metric report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*ALIGNMENT_SUMMARY_READ_2_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","PCT_PF_READS_ALIGNED_R2","PF_HQ_ALIGNED_READS_R2",\
"PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2",\
"PCT_READS_ALIGNED_IN_PAIRS_R2","PCT_ADAPTER_R2"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS_HEADER.TXT

# Sorting concatenated pair alignment summary metric report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*ALIGNMENT_SUMMARY_READ_PAIR_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","TOTAL_READS","RAW_GIGS","PCT_PF_READS_ALIGNED_PAIR",\
"PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR",\
"PCT_READS_ALIGNED_IN_PAIRS_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS_HEADER.TXT

# Sorting concatenated pair mark duplicates report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*MARK_DUPLICATES_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","UNMAPPED_READS",\
"READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS_HEADER.TXT

# Sorting concatenated pair hyb selection report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*HYB_SELECTION.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY",\
"PCT_PF_UQ_READS_ALIGNED","PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE",\
"ZERO_CVG_TARGETS_PCT","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET",\
"PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X",\
"AT_DROPOUT","GC_DROPOUT","HET_SNP_SENSITIVITY","HET_SNP_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION_HEADER.TXT

# Sorting concatenated bait bias report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*BAIT_BIAS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","Cref_Q","Gref_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS_HEADER.TXT

# Sorting concatenated pre adapter report by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*PRE_ADAPTER.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","DEAMINATION_Q","OxoG_Q"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER_HEADER.TXT

# Sorting concatenated ON BAIT SNV metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*BAIT_SNV_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_SNV_ON_BAIT","PERCENT_SNV_ON_BAIT_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/BAIT_SNV_METRICS_HEADER.TXT

# Sorting concatenated ON TARGET SNV metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TARGET_SNV_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_SNV_ON_TARGET","PERCENT_SNV_ON_TARGET_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TARGET_SNV_METRICS_HEADER.TXT

# Sorting concatenated TITV FOR ALL CODING SNVS metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TITV_ALL.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","ALL_TI_TV_COUNT","ALL_TI_TV_RATIO"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TITV_ALL_HEADER.TXT

# Sorting concatenated TITV FOR KNOWN CODING SNVS metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TITV_KNOWN.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","KNOWN_TI_TV_COUNT","KNOWN_TI_TV_RATIO"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TITV_KNOWN_HEADER.TXT

# Sorting concatenated TITV FOR NOVEL CODING SNVS metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TITV_NOVEL.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","NOVEL_TI_TV_COUNT","NOVEL_TI_TV_RATIO"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TITV_NOVEL_HEADER.TXT

# Sorting concatenated ON BAIT INDEL metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*BAIT_INDEL_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_ALL_INDEL_BAIT","ALL_INDEL_BAIT_PCT_SNP138","COUNT_BIALLELIC_INDEL_BAIT","BIALLELIC_INDEL_BAIT_PCT_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/BAIT_INDEL_METRICS_HEADER.TXT
# Sorting concatenated ON BAIT INDEL metrics by SM TAG and adding headers

# Sorting concatenated ON TARGET INDEL metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TARGET_INDEL_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_ALL_INDEL_TARGET","ALL_INDEL_TARGET_PCT_SNP138","COUNT_BIALLELIC_INDEL_TARGET","BIALLELIC_INDEL_TARGET_PCT_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TARGET_INDEL_METRICS_HEADER.TXT

# Sorting concatenated ON BAIT MIXED metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*BAIT_MIXED_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_MIXED_ON_BAIT","PERCENT_MIXED_ON_BAIT_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/BAIT_MIXED_METRICS_HEADER.TXT

# Sorting concatenated ON TARGET MIXED metrics by SM TAG and adding headers

cat $CORE_PATH/$PROJECT/TEMP/*TARGET_MIXED_METRICS.TXT \
| sort -k 1,1 \
| uniq \
| awk 'BEGIN {print "SM_TAG","COUNT_MIXED_ON_TARGET","PERCENT_MIXED_ON_TARGET_SNP138"} {print $0}' \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/TARGET_MIXED_METRICS_HEADER.TXT

######################################################################################################

# Joining all of the files together to make a QC report
# TO DO CONCATENATE WITH THE CONTROLS

TIMESTAMP=`date '+%F.%H-%M-%S'`s

join -i -1 2 -2 1 $CORE_PATH/$PROJECT/TEMP/SAMPLE_META_HEADER.txt $CORE_PATH/$PROJECT/TEMP/GENDER_CHECK_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/VERFIY_BAM_ID_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/INSERT_SIZE_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_1_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_2_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/ALIGNMENT_SUMMARY_READ_PAIR_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/MARK_DUPLICATES_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/HYB_SELECTION_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/BAIT_BIAS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/PRE_ADAPTER_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/BAIT_SNV_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TARGET_SNV_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TITV_ALL_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TITV_KNOWN_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TITV_NOVEL_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/BAIT_INDEL_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TARGET_INDEL_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/BAIT_MIXED_METRICS_HEADER.TXT \
| join -i -j 1 /dev/stdin $CORE_PATH/$PROJECT/TEMP/TARGET_MIXED_METRICS_HEADER.TXT \
| sed 's/ /,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".QC_REPORT."$TIMESTAMP".csv"

# Don't forget to add MD5s

# Concatenate all aneuploidy reports together
# TO DO CONCATENATE WITH THE CONTROLS...

( cat $CORE_PATH/$PROJECT/*/*/REPORTS/ANEUPLOIDY_CHECK/*.chrom_count_report.txt | grep "^SM_TAG" | uniq ; \
cat $CORE_PATH/$PROJECT/*/*/REPORTS/ANEUPLOIDY_CHECK/*.chrom_count_report.txt | grep -v "SM_TAG" ) \
| sed 's/\t/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".ANEUPLOIDY_CHECK."$TIMESTAMP".csv"

# Concatenate all per chromosome verifybamID reports together
# TO DO CONCATENATE WITH THE CONTROLS...

( cat $CORE_PATH/$PROJECT/*/*/REPORTS/VERIFYBAMID_CHR/*.VERIFYBAMID.PER_CHR.txt | grep "^#" | uniq ; \
cat $CORE_PATH/$PROJECT/*/*/REPORTS/VERIFYBAMID_CHR/*.VERIFYBAMID.PER_CHR.txt | grep -v "^#" ) \
| sed 's/\t/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".PER_CHR_VERIFYBAMID."$TIMESTAMP".csv"

# Summarize Wall Clock times

sed 's/,/\t/g' $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv" \
| sort -k 1,1 -k 2,2 -k 3,3 \
| awk 'BEGIN {OFS="\t"} {print $0,($6-$5),($6-$5)/60,($6-$5)/3600}' \
| $DATAMASH/datamash -s -g 1,2 max 7 max 8 max 9 | tee $CORE_PATH/$PROJECT/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
| $DATAMASH/datamash -g 1 sum 3 sum 4 sum 5 \
| awk 'BEGIN {print "SAMPLE_PROJECT","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
| sed -r 's/[[:space:]]+/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.BY_SAMPLE.csv"

sed 's/\t/,/g' $CORE_PATH/$PROJECT/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
| awk 'BEGIN {print "SAMPLE_PROJECT","TASK_GROUP","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
| sed -r 's/[[:space:]]+/,/g' \
>| $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.BY_SAMPLE_GROUP.csv"

md5sum $CORE_PATH/$PROJECT/LOGS/*log \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/*/LOGS/*log \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/*/*/LOGS/*log \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

echo Project finished at `date` >> $CORE_PATH/$PROJECT/REPORTS/PROJECT_START_END_TIMESTAMP.txt

# chmod 770 -R $CORE_PATH/$PROJECT
