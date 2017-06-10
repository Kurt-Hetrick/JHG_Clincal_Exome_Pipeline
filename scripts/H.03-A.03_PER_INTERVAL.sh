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

CORE_PATH=$1
BEDTOOLS_DIR=$2
CODING_BED=$3

PROJECT=$4
FAMILY=$5
SM_TAG=$6

# PAD THE ANNOTATED CODING BED FILE
# The input bed file could be a variable name based on the padding length

awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10,$4,$5,$6,$7}' $CODING_BED \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOTATED_PADDED_CODING_2.bed"

START_PER_INTERVAL=`date '+%s'`

awk 'NR>1' $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_summary.csv" \
| sed 's/:/\t/g; s/,/\t/g' \
| awk 'BEGIN {OFS="\t"} $2!~"-" {print $1,$2"-"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14} $2~"-" {print $0}' \
| awk 'BEGIN {OFS="\t"} {gsub (/-/,"\t"); print $1,$2,$3,$4,$5,$8,$9,$10,$13,$14,$15}' \
| $BEDTOOLS_DIR/bedtools intersect \
-a $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOTATED_PADDED_CODING_2.bed" \
-b - \
-wb \
| awk '{OFS="\t"} BEGIN {print "#CHROM" "\t" "START" "\t" "END" "\t" "GENE" "\t" "TRANSCRIPT" "\t" "EXON" "\t" "STRAND" "\t" "TOTAL_CVG" "\t" "AVG_CVG" "\t" \
"Q1_CVG" "\t" "MEDIAN_CVG" "\t" "Q3_CVG" "\t" "PCT_20x" "\t" "PCT_30x" "\t" "PCT_50x"} \
{print $1,$2,$3,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18}' \
>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.CODING_EXON_SUMMARY.txt"

END_PER_INTERVAL=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,REFSEQ_PER_INTERVAL,"$HOSTNAME","$START_PER_INTERVAL","$END_PER_INTERVAL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

md5sum $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.CODING_EXON_SUMMARY.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT.CIDR.Analysis.MD5.txt
