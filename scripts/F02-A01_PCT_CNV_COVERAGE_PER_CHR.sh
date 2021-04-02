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

	CNV_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	FAMILY=$4
	SM_TAG=$5

	CNV_CALL_PCT_BED=$6 # example command to create bed file
		# sed '1d' RefGene.Coding.Clean.Format.19Feb2018_gene_exomeDepth_noY.bed |sort -k1,1 -k2,2n -k3,3n -u|bedtools merge -i - >RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region.bed
		# BEDFILE=/mnt/research/statgen/DDL/CNV_pipeline_DDL/data/RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region2.bed
		# I think I know what this is...I also think it should be created from the original bed file.

	BED=(7821961 5837497 4499324 3091846 3506196 3878112 3658664 2566853 3146790 3101166 4523961 4255702 1432634 2414646 2787139 3189648 4417001 1220227 4581852 1857126 789143 1626513)

		# PENG'S COMMENT: for X it's 1297674, not inside BED. this is for the lengths of the chromosome.
		# KNH: this is the sum of the total length of intervals per chromosome in $CNV_CALL_PCT_BED.
		# KNH: I don't like this being hard-coded...
		# I would have this array be generated dynamically...
		# could do it here, but I think I would prefer
		# maybe put something in FIX_BED_FILES.sh to generate this array for autosomes
		# and have it dynamically generated here for chr X and chr Y
		# I think I'm leaning towards doing it here to create variables which then get passed.

# calculate the total length of cnv calls per autosome and divide by total length possible

for i in `seq 22`
do
	# take the exomeDepth output bam file
	awk -v chromosome="$i" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".exomeDepth.bed" \
	| sort -k 2,2n \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		merge -i - \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		intersect \
		-a - \
		-b $CNV_CALL_PCT_BED \
	| awk -v len=${BED[$i-1]} -v ch=$i '{SUM += $3-$2} END {printf "chr" ch " " "%.1f\n", 100*SUM/len}'
done \
>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".call.pct.by.chr.txt"

# calculate the total length of cnv calls for chrX and divide by total length possible

	awk -v chromosome="X" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".exomeDepth.bed" \
	| sort -k 2,2n \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		merge \
		-i - \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		intersect \
		-a - \
		-b $CNV_CALL_PCT_BED \
	| awk -v len="2760475" '{SUM += $3-$2} END {printf "chrX " "%.1f\n", 100*SUM/len}' \
	>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".call.pct.by.chr.txt"

# calculate the total length of cnv calls for chrX and divide by total length possible

	awk -v chromosome="Y" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".exomeDepth.bed" \
	| sort -k 2,2n \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		merge -i - \
	| singularity exec $CNV_CONTAINER \
		bedtools \
		intersect \
		-a - \
		-b $CNV_CALL_PCT_BED \
	| awk -v len="158436" '{SUM += $3-$2} END {printf "chrY " "%.1f\n", 100*SUM/len}' \
	>> $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT/$SM_TAG".call.pct.by.chr.txt"
