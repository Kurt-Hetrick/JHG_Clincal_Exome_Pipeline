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
	CODING_BED=$6
		CODING_BED_NAME=$(basename $CODING_BED .bed)
		CODING_MD5=$(md5sum $CODING_BED | cut -c 1-7)
	PADDING_LENGTH=$7

# Annotate the sample_interval_summary file from DEPTH_OF_COVERAGE.csv with what chromosome arm that intervals falls into
# remove any snps if present (they don't have a (-) in the first field, ignore the header
# split the first field on semicolon.
	# Now have two elements
		# First is CHROMOSOME
		# Second is the intervals. start(-)stop
		# Split second element on -
			# first element is START
			# second element is END
# print a bed FILE with tab-delimiter
	# CHROMOSOME
	# START-1 (to convert start to 0-based)
	# END
	# total_coverage (2nd field of file when using comma delimters)
# intersect said "bed file" with $SM_TAG".CHROM_ARM.bed" so that each record is annotated with what chromsome arm (p or q) it falls on
	# $SM_TAG".CHROM_ARM.bed" is created in "FIX_BED_FILES.sh" from the ucsc cytoband bed file
	# this is the part with bedtools intersect
	# also generates the number of bases that overlap (which is good enough so say how long the interval is)
# Reannotate CHROMOSOME X TARGETS TO X.PAR IF THEY FALL IN THE PSEUDOAUTOSOMAL REGIONS.
# If the CHROMOSOME is X and START is <=26995290 then change CHROMOSOME from X to X.PAR
	# print SM_TAG, X.PAR, ARM, total_coverage, interval length
# If the CHROMOSOME is X and START >=154931044 then change CHROMOSOME from X to X.PAR
	# print SM_TAG, X.PAR, ARM, total_coverage, interval length
# if neither condition above is met then just print the original CHROMOSOME
	# print SM_TAG, X.PAR, ARM, total_coverage, interval length
# group by SM_TAG,CHROMOSOME,ARM and sum up total_coverage and interval lengths
# group by SM_TAG,CHROMOSOME and sum up total coverage and interval lengths
	# now we have statistics (total number of sequenced bases for all intervals per whole chromosome as well as per chromosome arm)
	# and the length of bases attempted to be captured by whole chromosome and per chromosome arm

# Below does calculate for X PAR, but I don't think that it removes PAR from the X calculation...

	awk 'BEGIN {FS=","};{OFS="\t"} $1~"-" \
		{split($1,CHROM,":"); split(CHROM[2],POS,"-"); \
		print CHROM[1],POS[1]-1,POS[2],$2}' \
	$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/$SM_TAG"_"$CODING_BED_NAME"-"$CODING_MD5"-"$PADDING_LENGTH"-BP-PAD.sample_interval_summary.csv" \
		| singularity exec $ALIGNMENT_CONTAINER bedtools \
			intersect \
			-wo \
			-a - \
			-b $CORE_PATH/$PROJECT/TEMP/$SM_TAG".CHROM_ARM.bed" \
		| awk 'BEGIN {OFS="\t"} {if ($1=="X"&&$2<=2699520) print "'$SM_TAG'","X.PAR",$8,$4,$9 ; \
			else if ($1=="X"&&$2>=154931044) print "'$SM_TAG'","X.PAR",$8,$4,$9 ; \
			else print "'$SM_TAG'",$1,$8,$4,$9}' \
		| singularity exec $ALIGNMENT_CONTAINER \
			datamash \
			-s \
			-g 1,2,3 \
			sum 4 \
			sum 5 \
		| tee $CORE_PATH/$PROJECT/TEMP/$SM_TAG".depth_per_chr_arm.txt" \
		| singularity exec $ALIGNMENT_CONTAINER datamash \
			-g 1,2 \
			sum 4 \
			sum 5 \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,"whole",$3,$4}' \
		| tee -a $CORE_PATH/$PROJECT/TEMP/$SM_TAG".depth_per_chr_arm.txt"

# calculate the mean autosomal (chr 1-22) read depth to normalized all of the depth for each chromosome and chromsome arm

	AUTOSOMAL_MEAN_DEPTH=$(awk '$3=="whole"&&$2~/[0-9]/ {print $0}' \
		$CORE_PATH/$PROJECT/TEMP/$SM_TAG".depth_per_chr_arm.txt" \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				-s \
				sum 4 \
				sum 5 \
			| awk '{print $1/$2}')

# take the total coverages by chrom and chromosome arm
# calcuate the mean depth for each sample,chr,arm combination
# normalize by the AUTOSOMAL MEAN DEPTH

	awk 'BEGIN {print "SM_TAG","CHROM","ARM","TOTAL_COVERAGE","TOTAL_TARGETS","MEAN_DEPTH","NORM_DEPTH"} \
	{print $1,$2,$3,$4,$5,$4/$5,$4/$5/"'$AUTOSOMAL_MEAN_DEPTH'"}' \
	$CORE_PATH/$PROJECT/TEMP/$SM_TAG".depth_per_chr_arm.txt" \
		| awk '$2!="21"||$3!="p" {print $0}' \
		| awk '$2!="X.PAR"||$3!="p" {print $0}' \
		| awk '$2!="X.PAR"||$3!="q" {print $0}' \
		| sed 's/ /\t/g' \
	>| $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK/$SM_TAG".chrom_count_report.txt"
