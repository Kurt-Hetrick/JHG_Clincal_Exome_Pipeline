#!/usr/bin/env bash

# INPUT VARIABLES

	SAMPLE_SHEET=$1
	PED_FILE=$2
	PADDING_LENGTH=$3 # optional. if no 3rd argument present then the default is 10
	# THIS PAD IS FOR SLICING

		if [[ ! $PADDING_LENGTH ]]
			then
			PADDING_LENGTH="10"
		fi

	QUEUE_LIST=$4 # optional. if no 4th argument present then the default is cgc.q
		# if you want to set this then you need to set the 3rd argument as well (even to the default)

		if [[ ! $QUEUE_LIST ]]
			then
			QUEUE_LIST="cgc.q"
		fi

	PRIORITY=$5 # optional. if no 5th argument present then the default is -15.
		# if you want to set this then you need to set the 3rd and 4th argument as well (even to the default)

			if [[ ! $PRIORITY ]]
				then
				PRIORITY="-15"
			fi

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

	SUBMITTER_SCRIPT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

	SCRIPT_DIR="$SUBMITTER_SCRIPT_PATH/scripts"

##################
# CORE VARIABLES #
##################

	# GVCF PAD. CURRENTLY KEEPING THIS AS A STATIC VARIABLE

		GVCF_PAD="250"

	## This will always put the current working directory in front of any directory for PATH
	## added /bin for RHEL6

		export PATH=".:$PATH:/bin"

	# where the input/output sequencing data will be located.

		CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"

	# Directory where NovaSeqa runs are located.

		NOVASEQ_REPO="/mnt/instrument_files/novaseq"

	# this is just for note taking where the control data set will reside. not used in this pipeline.

		CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Resources/TWIST_CONTROL_SET1.200601_PIPELINE_2_0_0"

	# used for tracking in the read group header of the cram file

		PIPELINE_VERSION=`git --git-dir=$SCRIPT_DIR/../.git --work-tree=$SCRIPT_DIR/.. log --pretty=format:'%h' -n 1`

	# load gcc for programs like verifyBamID
	## this will get pushed out to all of the compute nodes since I specify env var to pushed out with qsub

		module load gcc/7.2.0

	# explicitly setting this b/c not everybody has had the $HOME directory transferred and I'm not going to through
	# and figure out who does and does not have this set correctly

		umask 0007

	# SUBMIT TIMESTAMP

		SUBMIT_STAMP=`date '+%s'`

	# SUBMITTER_ID

		SUBMITTER_ID=`whoami`

	# bind the host file system /mnt to the singularity container. in case I use it in the submitter.

		export SINGULARITY_BINDPATH="/mnt:/mnt"

	# QSUB ARGUMENTS LIST
		# set shell on compute node
		# start in current working directory
		# transfer submit node env to compute node
		# set SINGULARITY BINDPATH
		# set queues to submit to
		# set priority
		# combine stdout and stderr logging to same output file

			QSUB_ARGS="-S /bin/bash" \
				QSUB_ARGS=$QSUB_ARGS" -cwd" \
				QSUB_ARGS=$QSUB_ARGS" -V" \
				QSUB_ARGS=$QSUB_ARGS" -v SINGULARITY_BINDPATH=/mnt:/mnt" \
				QSUB_ARGS=$QSUB_ARGS" -q $QUEUE_LIST" \
				QSUB_ARGS=$QSUB_ARGS" -p $PRIORITY" \
				QSUB_ARGS=$QSUB_ARGS" -j y"

#####################
# PIPELINE PROGRAMS #
#####################

	ALIGNMENT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/images/ddl_ce_control_align-0.0.3.simg"
	# contains the following software and is on Ubuntu 16.04.5 LTS
		# gatk 4.0.11.0 (base image). also contains the following.
			# Python 3.6.2 :: Continuum Analytics, Inc.
				# samtools 0.1.19
				# bcftools 0.1.19
				# bedtools v2.25.0
				# bgzip 1.2.1
				# tabix 1.2.1
				# samtools, bcftools, bgzip and tabix will be replaced with newer versions.
				# R 3.2.5
					# dependencies = c("gplots","digest", "gtable", "MASS", "plyr", "reshape2", "scales", "tibble", "lazyeval")    # for ggplot2
					# getopt_1.20.0.tar.gz
					# optparse_1.3.2.tar.gz
					# data.table_1.10.4-2.tar.gz
					# gsalib_2.1.tar.gz
					# ggplot2_2.2.1.tar.gz
				# openjdk version "1.8.0_181"
				# /gatk/gatk.jar -> /gatk/gatk-package-4.0.11.0-local.jar
		# added
			# picard.jar 2.17.0 (as /gatk/picard.jar)
			# samblaster-v.0.1.24
			# sambamba-0.6.8
			# bwa-0.7.15
			# datamash-1.6
			# verifyBamID v1.1.3
			# samtools 1.10
			# bgzip 1.10
			# tabix 1.10
			# bcftools 1.10.2

	GATK_3_7_0_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/images/gatk3-3.7-0.simg"
	# singularity pull docker://broadinstitute/gatk3:3.7-0
	# used for generating the depth of coverage reports.
		# comes with R 3.1.1 with appropriate packages needed to create gatk pdf output
		# also comes with some version of java 1.8
		# jar file is /usr/GenomeAnalysisTK.jar

	# PIPELINE PROGRAMS TO BE IMPLEMENTED
	JAVA_1_6="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jre1.6.0_25/bin"
	SAMTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/samtools-0.1.18"
	VCFTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/vcftools_0.1.12b/bin"
	PLINK2_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/PLINK2"
	KING_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/KING/Linux-king19"
	CIDRSEQSUITE_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/CIDRSeqSuiteSoftware_Version_4_0/"
	ANNOVAR_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/ANNOVAR/2013_09_11"

##################
# PIPELINE FILES #
##################

	GENE_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/RefSeqGene.GRCh37.rCRS.MT.bed"
		# md5 dec069c279625cfb110c2e4c5480e036
	VERIFY_VCF="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
	CODING_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINES/TWIST/JHGenomics_CGC_Clinical_Exome_Control_Set/GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_annotated.bed"
	CYTOBAND_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/GRCh37.Cytobands.bed"
	HAPMAP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/hapmap_3.3.b37.vcf"
	OMNI_1KG="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
	HI_CONF_1KG_PHASE1_SNP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
	MILLS_1KG_GOLD_INDEL="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
	PHASE3_1KG_AUTOSOMES="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
	DBSNP_129="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"
	CONTROL_PED_FILE="$CONTROL_REPO/TWIST_CONTROL_SET1.200601.ped"

#################################
##### MAKE A DIRECTORY TREE #####
#################################

##### CREATING A DIRECTORY IN USER'S HOME TO MERGE THE SAMPLE MANIFEST WITH THE PED FILE

	mkdir -p ~/CGC_PIPELINE_TEMP

	MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
	PED_PREFIX=`basename $PED_FILE .ped`

	FORMAT_MANIFEST ()
	{
		awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
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
		| join -1 8 -2 2 -e '-'  -t $'\t' \
		-o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
		~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
		>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
	}

	CREATE_SAMPLE_ARRAY ()
	{
		SAMPLE_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" \
			{split($19,INDEL,";"); \
			print $1,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$20,$21,$22,$23,$24}' \
				~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
				| sort \
				| uniq`)

		#  1  Project=the Seq Proj folder name

			PROJECT=${SAMPLE_ARRAY[0]}

		################################################################################
		# 2 SKIP : FCID=flowcell that sample read group was performed on ###############
		# 3 SKIP : Lane=lane of flowcell that sample read group was performed on] ######
		# 4 SKIP : Index=sample barcode ################################################
		# 5 SKIP : Platform=type of sequencing chemistry matching SAM specification ####
		# 6 SKIP : Library_Name=library group of the sample read group #################
		# 7 SKIP : Date=should be the run set up date to match the seq run folder name #
		################################################################################

		#  8  SM_Tag=sample ID

			SM_TAG=${SAMPLE_ARRAY[1]}
				SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # If there is an @ in the qsub or holdId name it breaks

		#  9  Center=the center/funding mechanism

			CENTER=${SAMPLE_ARRAY[2]}

		# 10  Description=Generally we use to denote the sequencer setting (e.g. rapid run)
		# “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.

			SEQUENCER_MODEL=${SAMPLE_ARRAY[3]}

		#########################
		# 11  SKIP : Seq_Exp_ID #
		#########################

		# 12  Genome_Ref=the reference genome used in the analysis pipeline

			REF_GENOME=${SAMPLE_ARRAY[4]}

		#####################################
		# 13  Operator: SKIP ################
		# 14  Extra_VCF_Filter_Params: SKIP #
		#####################################

		# 15  TS_TV_BED_File=where ucsc coding exons overlap with bait and target bed files

			TITV_BED=${SAMPLE_ARRAY[5]}

		# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
		# Used for limited where to run base quality score recalibration on where to create gvcf files.

			BAIT_BED=${SAMPLE_ARRAY[6]}

		# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.

			TARGET_BED=${SAMPLE_ARRAY[7]}

		# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.

			DBSNP=${SAMPLE_ARRAY[8]}

		# 19  KNOWN_INDEL_FILES=used for BQSR masking, sensitivity in local realignment.

			KNOWN_INDEL_1=${SAMPLE_ARRAY[9]}
			KNOWN_INDEL_2=${SAMPLE_ARRAY[10]}

		# 20 family that sample belongs to

			FAMILY=${SAMPLE_ARRAY[11]}

		# 21 MOM

			MOM=${SAMPLE_ARRAY[12]}

		# 22 DAD

			DAD=${SAMPLE_ARRAY[13]}

		# 23 GENDER

			GENDER=${SAMPLE_ARRAY[14]}

		# 24 PHENOTYPE

			PHENOTYPE=${SAMPLE_ARRAY[15]}
	}

# PROJECT DIRECTORY TREE CREATOR

	MAKE_PROJ_DIR_TREE ()
	{
		mkdir -p $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/HC_CRAM \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/INDEL/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/SNV/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MIXED/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_CHR} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/CONCORDANCE \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/{TARGET,REFSEQ_CODING_PLUS_10bp} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK \
		$CORE_PATH/$PROJECT/$FAMILY/{LOGS,VCF,RELATEDNESS,PCA} \
		$CORE_PATH/$PROJECT/TEMP/$SM_TAG_ANNOVAR \
		$CORE_PATH/$PROJECT/TEMP/{VCF_PREP,PLINK,KING} \
		$CORE_PATH/$PROJECT/{TEMP,FASTQ,REPORTS,LOGS,COMMAND_LINES}
	}

	SETUP_PROJECT ()
	{
		FORMAT_MANIFEST
		MERGE_PED_MANIFEST
		CREATE_SAMPLE_ARRAY
		MAKE_PROJ_DIR_TREE
		echo Project started at `date` >| $CORE_PATH/$PROJECT/REPORTS/PROJECT_START_END_TIMESTAMP.txt
	}

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		SETUP_PROJECT
done

# MAKE AN ARRAY FOR EACH SAMPLE
	## SAMPLE_INFO_ARRAY[0] = PROJECT
	## SAMPLE_INFO_ARRAY[1] = FAMILY
	## SAMPLE_INFO_ARRAY[2] = SM_TAG
		## SAMPLE = SM_TAG
	## SAMPLE_INFO_ARRAY[3] = BAIT BED FILE
	## SAMPLE_INFO_ARRAY[4] = TARGET_BED_FILE

############################################################

################################
##### CRAM FILE GENERATION #####
###############################################################################################
##### NOTE: THE CRAM FILE IS THE END PRODUCT BUT THE BAM FILE IS USED FOR OTHER PROCESSES #####
##### SOME PROGRAMS CAN'T TAKE IN CRAM AS AN INPUT ############################################
###############################################################################################

	########################################################################################
	# create an array at the platform level so that bwa mem can add metadata to the header #
	########################################################################################

		CREATE_PLATFORM_UNIT_ARRAY ()
		{
			PLATFORM_UNIT_ARRAY=(`awk 1 ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS="\t"} $8$2$3$4=="'$PLATFORM_UNIT'" {split($19,INDEL,";"); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$$20,$21,$22,$23,$24}' \
			| sort \
			| uniq`)

				#  1  Project=the Seq Proj folder name

					PROJECT=${PLATFORM_UNIT_ARRAY[0]}

				#  2  FCID=flowcell that sample read group was performed on

					FCID=${PLATFORM_UNIT_ARRAY[1]}

				#  3  Lane=lane of flowcell that sample read group was performed on

					LANE=${PLATFORM_UNIT_ARRAY[2]}

				#  4  Index=sample barcode

					INDEX=${PLATFORM_UNIT_ARRAY[3]}

				#  5  Platform=type of sequencing chemistry matching SAM specification

					PLATFORM=${PLATFORM_UNIT_ARRAY[4]}

				#  6  Library_Name=library group of the sample read group,
					# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not

					LIBRARY=${PLATFORM_UNIT_ARRAY[5]}

				#  7  Date=should be the run set up date, but doesn't have to be

					RUN_DATE=${PLATFORM_UNIT_ARRAY[6]}

				#  8  SM_Tag=sample ID

					SM_TAG=${PLATFORM_UNIT_ARRAY[7]}

						# sge sm tag. If there is an @ in the qsub or holdId name it breaks

							SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g')

				#  9  Center=the center/funding mechanism

					CENTER=${PLATFORM_UNIT_ARRAY[8]}

				# 10  Description=Generally we use to denote the sequencer setting (e.g. rapid run)
				# “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.

					SEQUENCER_MODEL=${PLATFORM_UNIT_ARRAY[9]}

				########################
				# 11  Seq_Exp_ID: SKIP #
				########################

				# 12  Genome_Ref=the reference genome used in the analysis pipeline

					REF_GENOME=${PLATFORM_UNIT_ARRAY[10]}

				#####################################
				# 13  Operator: SKIP ################
				# 14  Extra_VCF_Filter_Params: SKIP #
				#####################################

				# 15  TS_TV_BED_File=refseq (select) cds plus other odds and ends (.e.g. missing omim))

					TITV_BED=${PLATFORM_UNIT_ARRAY[11]}

				# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
				# Used for limited where to run base quality score recalibration on where to create gvcf files.

					BAIT_BED=${PLATFORM_UNIT_ARRAY[12]}

				# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.

					TARGET_BED=${PLATFORM_UNIT_ARRAY[13]}

				# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.

					DBSNP=${PLATFORM_UNIT_ARRAY[14]}

				# 19  KNOWN_INDEL_FILES=used for BQSR masking

					KNOWN_INDEL_1=${PLATFORM_UNIT_ARRAY[15]}
					KNOWN_INDEL_2=${PLATFORM_UNIT_ARRAY[16]}

				# 20 FAMILY

					FAMILY=${PLATFORM_UNIT_ARRAY[17]}

				# 21 MOM

					MOM=${PLATFORM_UNIT_ARRAY[17]}

				# 22 DAD

					DAD=${PLATFORM_UNIT_ARRAY[17]}

				# 23 GENDER

					GENDER=${PLATFORM_UNIT_ARRAY[17]}

				# 24 PHENOTYPE

					PHENOTYPE=${PLATFORM_UNIT_ARRAY[17]}
		}

	########################################################################
	### Use bwa mem to do the alignments; ##################################
	### pipe to samblaster to add mate tags; ###############################
	### pipe to picard's AddOrReplaceReadGroups to handle the bam header ###
	########################################################################

		RUN_BWA ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N A.01-BWA"_"$SGE_SM_TAG"_"$FCID"_"$LANE"_"$INDEX \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"_"$FCID"_"$LANE"_"$INDEX"-BWA.log" \
			$SCRIPT_DIR/A.01_BWA.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FCID \
				$LANE \
				$INDEX \
				$PLATFORM \
				$LIBRARY \
				$RUN_DATE \
				$SM_TAG \
				$CENTER \
				$SEQUENCER_MODEL \
				$REF_GENOME \
				$PIPELINE_VERSION \
				$BAIT_BED \
				$TARGET_BED \
				$TITV_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP \
				$NOVASEQ_REPO
		}

	for PLATFORM_UNIT in $(awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $8$2$3$4}' \
			| sort \
			| uniq );
		do
			CREATE_PLATFORM_UNIT_ARRAY
			mkdir -p $CORE_PATH/$PROJECT/LOGS/$SM_TAG
			RUN_BWA
			echo sleep 0.1s
	done

	#########################################################################################
	# Merge files and mark duplicates using picard duplictes with queryname sorting #########
	# do coordinate sorting with sambamba ###################################################
	#########################################################################################
	# I am setting the heap space and garbage collector threads for picard now now ##########
	# doing this does drastically decrease the load average ( the gc thread specification ) #
	#########################################################################################
	# create a hold job id qsub command line based on the number of #########################
	# submit merging the bam files created by bwa mem above #################################
	# only launch when every lane for a sample is done being processed by bwa mem ###########
	# I want to clean this up eventually and get away from using awk to print the qsub line #
	#########################################################################################

	# What is being pulled out of the merged sample sheet and ped file table.
		# 1. PROJECT
		# 2. FAMILY
		# 3. SM_TAG
		# 4. FCID_LANE_INDEX
		# 5. FCID_LANE_INDEX.bam
		# 6. SM_TAG
		# 7. DESCRIPTION (INSTRUMENT MODEL)

			awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam",$8,$10}' \
			~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
			| awk 'BEGIN {OFS="\t"} {sub(/@/,"_",$6)} {print $1,$2,$3,$4,$5,$6,$7}' \
			| sort -k 1,1 -k 2,2 -k 3,3 -k 4,4 -k 7,7 \
			| uniq \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				-s \
				-g 1,2,3 \
				collapse 4 \
				collapse 5 \
				unique 6 \
				unique 7 \
			| awk 'BEGIN {FS="\t"} \
				gsub(/,/,",A.01_BWA_"$6"_",$4) \
				gsub(/,/,",INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/",$5) \
				{print "qsub",\
				"-S /bin/bash",\
				"-cwd",\
				"-V",\
				"-v SINGULARITY_BINDPATH=/mnt:/mnt",\
				"-q","'$QUEUE_LIST'",\
				"-p","'$PRIORITY'",\
				"-j y",\
				"-N","B.01-MARK_DUPLICATES_"$6"_"$1,\
				"-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1"-MARK_DUPLICATES.log",\
				"-hold_jid","A.01-BWA_"$6"_"$4, \
				"'$SCRIPT_DIR'""/B.01_MARK_DUPLICATES.sh",\
				"'$ALIGNMENT_CONTAINER'",\
				"'$CORE_PATH'",\
				$1,\
				$2,\
				$3,\
				$7,\
				"'$SAMPLE_SHEET'",\
				"'$SUBMIT_STAMP'",\
				"INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/"$5"\n""sleep 0.1s"}'

	###############################################
	# fix common formatting problems in bed files #
	# merge bait to target for gvcf creation, pad #
	# create picard style interval files ##########
	###############################################

		FIX_BED_FILES ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-FIX_BED_FILES.log" \
			-hold_jid B.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/C.01_FIX_BED.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$TARGET_BED \
				$BAIT_BED \
				$TITV_BED \
				$CYTOBAND_BED \
				$REF_GENOME \
				$PADDING_LENGTH \
				$GVCF_PAD
		}

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		CREATE_SAMPLE_ARRAY
		FIX_BED_FILES
		echo sleep 0.1s
		# PERFORM_BQSR
		# echo sleep 0.1s
		# APPLY_BQSR
		# echo sleep 0.1s
		# BAM_TO_CRAM
		# echo sleep 0.1s
		# INDEX_CRAM
		# echo sleep 0.1s
done

# # Run Base Quality Score Recalibration

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19,$18,$16}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
# print "qsub","-N","F.01_PERFORM_BQSR_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","E.01_INDEL_REALIGNER_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PERFORM_BQSR.log",\
# "'$SCRIPT_DIR'""/F.01_PERFORM_BQSR.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6,$7"\n""sleep 1s"}'

# # write Final Bam file

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","F.01_PERFORM_BQSR_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".FINAL_BAM.log",\
# "'$SCRIPT_DIR'""/G.01_FINAL_BAM.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# # SCATTER THE HAPLOTYPE CALLER GVCF CREATION USING THE WHERE THE BED INTERSECTS WITH {{1.22},{X,Y}}

# CREATE_SAMPLE_INFO_ARRAY_HC ()
# {
# SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_HAPLOTYPE_CALLER ()
# {
# echo \
# qsub \
# -N H.01_HAPLOTYPE_CALLER_${SAMPLE_INFO_ARRAY_HC[0]}_${SAMPLE_INFO_ARRAY_HC[5]}_chr$CHROMOSOME \
# -hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_HC[5]}_${SAMPLE_INFO_ARRAY_HC[0]} \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_chr$CHROMOSOME.log \
# $SCRIPT_DIR/H.01_HAPLOTYPE_CALLER_SCATTER.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]} \
# $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_HC
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_HAPLOTYPE_CALLER
# 		echo sleep 1s
# 		done
# 	done

# ################################################################

# # GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES INTO A SINGLE SAMPLE GVCF

# # GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES INTO A SINGLE SAMPLE GVCF

# # BUILD_HOLD_ID_PATH(){
# # 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# # 	do
# # 	HOLD_ID_PATH="-hold_jid "
# # 	for CHROMOSOME in {{1..22},{X,Y}};
# #  	do
# #  		HOLD_ID_PATH=$HOLD_ID_PATH"H.01_HAPLOTYPE_CALLER_"$PROJECT"_"${SAMPLE_INFO_ARRAY_HC[5]}"_chr"$CHROMOSOME","
# #  	done
# #  done
# # }

# BUILD_HOLD_ID_PATH(){
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"H.01_HAPLOTYPE_CALLER_"$PROJECT"_"$SAMPLE"_chr"$CHROMOSOME","
#  	done
#  done
# }

# # CREATE_SAMPLE_INFO_ARRAY_HC ()
# # {
# # SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# # ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# CREATE_SAMPLE_INFO_ARRAY_HC ()
# {
# SAMPLE_INFO_ARRAY_HC=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_HAPLOTYPE_CALLER_GATHER ()
# {
# echo \
# qsub \
# -N H.01-A.01_HAPLOTYPE_CALLER_GATHER_${SAMPLE_INFO_ARRAY_HC[0]}_$SAMPLE \
# ${HOLD_ID_PATH} \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_GATHER.log \
# $SCRIPT_DIR/H.01-A.01_HAPLOTYPE_CALLER_GATHER.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]}
# }

# # CALL_HAPLOTYPE_CALLER_GATHER ()
# # {
# # echo \
# # qsub \
# # -N H.01-A.01_HAPLOTYPE_CALLER_GATHER_${SAMPLE_INFO_ARRAY_HC[0]}_${SAMPLE_INFO_ARRAY_HC[5]} \
# # ${HOLD_ID_PATH} \
# # -o $CORE_PATH/${SAMPLE_INFO_ARRAY_HC[0]}/${SAMPLE_INFO_ARRAY_HC[1]}/${SAMPLE_INFO_ARRAY_HC[2]}/LOGS/${SAMPLE_INFO_ARRAY_HC[2]}_${SAMPLE_INFO_ARRAY_HC[0]}.HAPLOTYPE_CALLER_GATHER.log \
# # $SCRIPT_DIR/H.01-A.01_HAPLOTYPE_CALLER_GATHER.sh \
# # $JAVA_1_8 $GATK_DIR $CORE_PATH \
# # ${SAMPLE_INFO_ARRAY_HC[0]} ${SAMPLE_INFO_ARRAY_HC[1]} ${SAMPLE_INFO_ARRAY_HC[2]} ${SAMPLE_INFO_ARRAY_HC[3]}
# # }

# # for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8} $SAMPLE_SHEET | sort | uniq );
# #  do
# # 	BUILD_HOLD_ID_PATH
# # 	CREATE_SAMPLE_INFO_ARRAY_HC
# # 	CALL_HAPLOTYPE_CALLER_GATHER
# # 	echo sleep 1s
# #  done

# # for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]}' $SAMPLE_SHEET | sort | uniq );

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
#  do
# 	BUILD_HOLD_ID_PATH
# 	CREATE_SAMPLE_INFO_ARRAY_HC
# 	CALL_HAPLOTYPE_CALLER_GATHER
# 	echo sleep 1s
#  done

# # Run POST BQSR TABLE

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$19,$18}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
# print "qsub","-N","H.02_POST_BQSR_TABLE_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".POST_BQSR_TABLE.log",\
# "'$SCRIPT_DIR'""/H.02_POST_BQSR_TABLE.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,INDEL[1],INDEL[2],$6"\n""sleep 1s"}'

# # Run ANALYZE COVARIATES

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($5,INDEL,";"); split($3,smtag,"[@]"); \
# print "qsub","-N","H.02-A.01_ANALYZE_COVARIATES_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.02_POST_BQSR_TABLE_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANALYZE_COVARIATES.log",\
# "'$SCRIPT_DIR'""/H.02-A.01_ANALYZE_COVARIATES.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# # RUN DOC CODING PLUS 10 BP FLANKS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_CODING_10bpFLANKS.log",\
# "'$SCRIPT_DIR'""/H.03_DOC_CODING_10bpFLANKS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 1s"}'

# # RUN ANEUPLOIDY_CHECK AFTER CODING PLUS 10 BP FLANKS FINISHES

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.01_DOC_CHROM_DEPTH_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".ANEUPLOIDY_CHECK.log",\
# "'$SCRIPT_DIR'""/H.03-A.01_CHROM_DEPTH.sh",\
# "'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2,$3"\n""sleep 1s"}'

# # RUN FORMATTING PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE.log",\
# "'$SCRIPT_DIR'""/H.03-A.02_PER_BASE.sh",\
# "'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2,$3"\n""sleep 1s"}'

# # RUN FILTERING PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02_PER_BASE_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_FILTER.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.01_PER_BASE_FILTERED.sh",\
# "'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_BGZIP.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.02_PER_BASE_BGZIP.sh",\
# "'$CORE_PATH'","'$TABIX_DIR'",$1,$2,$3"\n""sleep 1s"}'

# # TABIX PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02-A.02-A.01_PER_BASE_TABIX_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_BASE_TABIX.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.02-A.01_PER_BASE_TABIX.sh",\
# "'$CORE_PATH'","'$TABIX_DIR'",$1,$2,$3"\n""sleep 1s"}'

# # RUN FORMATTING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03_DOC_CODING_10bpFLANKS_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_INTERVAL.log",\
# "'$SCRIPT_DIR'""/H.03-A.03_PER_INTERVAL.sh",\
# "'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2,$3"\n""sleep 1s"}'

# # RUN FILTERING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.03-A.03_PER_INTERVAL_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".PER_INTERVAL_FILTER.log",\
# "'$SCRIPT_DIR'""/H.03-A.03-A.01_PER_INTERVAL_FILTERED.sh",\
# "'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # RUN DOC TARGET BED

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$17}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.05_DOC_TARGET_BED_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".DOC_TARGET_BED.log",\
# "'$SCRIPT_DIR'""/H.05_DOC_TARGET_BED.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 1s"}'

# # RUN COLLECT MULTIPLE METRICS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$18,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.06_COLLECT_MULTIPLE_METRICS_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_MULTIPLE_METRICS.log",\
# "'$SCRIPT_DIR'""/H.06_COLLECT_MULTIPLE_METRICS.sh",\
# "'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5,$6"\n""sleep 1s"}'

# # RUN COLLECT HS METRICS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$16,$17}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.07_COLLECT_HS_METRICS_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".COLLECT_HS_METRICS.log",\
# "'$SCRIPT_DIR'""/H.07_COLLECT_HS_METRICS.sh",\
# "'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4"\n""sleep 1s"}'

# # RUN SELECT VERIFYBAM ID VCF

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.08_SELECT_VERIFYBAMID_VCF_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".SELECT_VERIFYBAMID_VCF.log",\
# "'$SCRIPT_DIR'""/H.08_SELECT_VERIFYBAMID_VCF.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VERIFY_VCF'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# # RUN VERIFYBAMID

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); \
# print "qsub","-N","H.08-A.01_VERIFYBAMID_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.08_SELECT_VERIFYBAMID_VCF_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$1".VERIFYBAMID.log",\
# "'$SCRIPT_DIR'""/H.08-A.01_VERIFYBAMID.sh",\
# "'$CORE_PATH'","'$VERIFY_DIR'",$1,$2,$3"\n""sleep 1s"}'

# ###################################################
# ### RUN VERIFYBAM ID PER CHROMOSOME - VITO ########
# ###################################################

# CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM ()
# {
# SAMPLE_INFO_ARRAY_VERIFY_BAM=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$20,$8,$12,$15,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_SELECT_VERIFY_BAM ()
# {
# echo \
# qsub \
# -N H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.SELECT_VERIFYBAMID_chr$CHROMOSOME.log \
# $SCRIPT_DIR/H.09_SELECT_VERIFYBAMID_VCF_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH $VERIFY_VCF \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[3]} \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[4]} $CHROMOSOME
# }

# CALL_VERIFYBAMID ()
# {
# echo \
# qsub \
# -N H.09-A.01_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -hold_jid H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.VERIFYBAMID_chr$CHROMOSOME.log \
# $SCRIPT_DIR/H.09-A.01_VERIFYBAMID_CHR.sh \
# $CORE_PATH $VERIFY_DIR \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} \
# $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
# 	for CHROMOSOME in {1..22}
# 		do
# 		CALL_SELECT_VERIFY_BAM
# 		echo sleep 1s
# 		CALL_VERIFYBAMID
# 		echo sleep 1s
# 	done
# done

# #####################################################
# ### JOIN THE PER CHROMOSOME VERIFYBAMID REPORTS #####
# #####################################################

# BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"H.09-A.01_VERIFYBAMID_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}"_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}"_"chr$CHROMOSOME","
#  	done
#  done
# }

#  CAT_VERIFYBAMID_CHR ()
#  {
# echo \
# qsub \
# -N H.09-A.01-A.01_JOIN_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[5]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
# $HOLD_ID_PATH \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]}/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.CAT_VERIFYBAMID_CHR.log \
# $SCRIPT_DIR/H.09-A.01-A.01_CAT_VERIFYBAMID_CHR.sh \
# $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[1]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}
#  }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
#  do
#  	CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
# 	BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR
# 	CAT_VERIFYBAMID_CHR
# 	echo sleep 1s
#  done

# #### JOINT CALLING AND VQSR #### ###VITO###

# ### CREATE A GVCF ".list" file for each sample

# CREATE_FAMILY_INFO_ARRAY ()
# {
# FAMILY_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {split($8,smtag,"[@]"); print $1,smtag[1]"_"smtag[2],$20,$12,$18}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CREATE_GVCF_LIST ()
# {
# awk 'BEGIN {OFS="/"} $20=="'$FAMILY'" {print "'$CORE_PATH'",$1,$20,$8,"GVCF",$8".g.vcf.gz"}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort \
# | uniq \
# >| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".gvcf.list"
# }

# CREATE_FAMILY_SAMPLE_LIST ()
# {
# awk '$20=="'$FAMILY'" {print $8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort \
# | uniq \
# >| $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/$FAMILY/$FAMILY".sample.list"
# }

# BUILD_HOLD_ID_PATH_GENOTYPE_GVCF ()
# {
# ##NEED FULL LIST OF SAMPLES IN FAMILY FOR THE INNER FOR LOOP.  THE ${FAMILY_INFO_ARRAY[1]} WON'T WORK CONSIDERING THAT IS SPECIFIC FOR ONE SAMPLE

# for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# do
# HOLD_ID_PATH="-hold_jid "
# for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
# 	do
# 		HOLD_ID_PATH=$HOLD_ID_PATH"H.01-A.01_HAPLOTYPE_CALLER_GATHER_"$PROJECT"_"${FAMILY_INFO_ARRAY[1]}","
# 	done
# done
# }

# CALL_GENOTYPE_GVCF ()
# {
# for CHROM in {{1..22},{X,Y}};
# do
# echo \
# qsub \
# -N I.01_GENOTYPE_GVCF_SCATTER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}_chr$CHROM \
# $HOLD_ID_PATH \
# -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.GENOTYPE_GVCF_chr$CHROM.log \
# $SCRIPT_DIR/I.01_GENOTYPE_GVCF_SCATTER.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]} ${FAMILY_INFO_ARRAY[4]} $CHROM $CONTROL_REPO
# echo sleep 1s
# done
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
# do
# 	 CREATE_FAMILY_INFO_ARRAY
# 	 CREATE_GVCF_LIST
#  CREATE_FAMILY_SAMPLE_LIST
# BUILD_HOLD_ID_PATH_GENOTYPE_GVCF
# CALL_GENOTYPE_GVCF
# done

# ########################################################################################
# ##### GATHER UP THE PER FAMILY PER CHROMOSOME GVCF FILES INTO A SINGLE FAMILY GVCF #####
# ########################################################################################

# BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHR in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"I.01_GENOTYPE_GVCF_SCATTER_"$FAMILY"_"$PROJECT"_chr"$CHR","
#  	done
#  done
# }


# CALL_GENOTYPE_GVCF_GATHER ()
# {
# echo \
# qsub \
# -N I.01-A.01_GENOTYPE_GVCF_GATHER_${FAMILY_INFO_ARRAY[0]}_${FAMILY_INFO_ARRAY[2]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.GENOTYPE_GVCF_GATHER.log \
#  $SCRIPT_DIR/I.01-A.01_GENOTYPE_GVCF_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }


# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
#  	# echo $FAMILY
# 	BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_GENOTYPE_GVCF_GATHER
# 	echo sleep 1s
#  done

# ##########################################################
# ################END VITO##################################
# ##########################################################

# #####################################################################################################
# ##### Run Variant Recalibrator for the SNP model, this is done in parallel with the INDEL model #####
# #####################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12,$18}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","J.01_VARIANT_RECALIBRATOR_SNP_"$2"_"$1,\
# "-hold_jid","I.01-A.01_GENOTYPE_GVCF_GATHER_"$1"_"$2,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_SNP.log",\
# "'$SCRIPT_DIR'""/J.01_VARIANT_RECALIBRATOR_SNP.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,"'$HAPMAP'","'$OMNI_1KG'","'$HI_CONF_1KG_PHASE1_SNP'""\n""sleep 1s"}'

# ### Run Variant Recalibrator for the INDEL model, this is done in parallel with the SNP model

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","J.02_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
# "-hold_jid","I.01-A.01_GENOTYPE_GVCF_GATHER_"$1"_"$2,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_RECALIBRATOR_INDEL.log",\
# "'$SCRIPT_DIR'""/J.02_VARIANT_RECALIBRATOR_INDEL.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,"'$MILLS_1KG_GOLD_INDEL'""\n""sleep 1s"}'

# ### Run Apply Recalbration with the SNP model to the VCF file

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","K.01_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
# "-hold_jid","J.01_VARIANT_RECALIBRATOR_SNP_"$2"_"$1",""J.02_VARIANT_RECALIBRATOR_INDEL_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_SNP.log",\
# "'$SCRIPT_DIR'""/K.01_APPLY_RECALIBRATION_SNP.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ### Run Apply Recalibration with the INDEL model to the VCF file.

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","L.01_APPLY_RECALIBRATION_INDEL_"$2"_"$1,\
# "-hold_jid","K.01_APPLY_RECALIBRATION_SNP_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".APPLY_RECALIBRATION_INDEL.log",\
# "'$SCRIPT_DIR'""/L.01_APPLY_RECALIBRATION_INDEL.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ################################################
# ##### SCATTER GATHER FOR ADDING ANNOTATION #####
# ################################################

# CREATE_FAMILY_INFO_ARRAY ()
# {
# FAMILY_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $1,$8,$20,$12,$18}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_VARIANT_ANNOTATOR ()
# {
# echo \
# qsub \
# -N P.01_VARIANT_ANNOTATOR_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}_$CHROMOSOME \
# -hold_jid L.01_APPLY_RECALIBRATION_INDEL_${FAMILY_INFO_ARRAY[2]}"_"${FAMILY_INFO_ARRAY[0]} \
# -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_ANNOTATOR_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01_VARIANT_ANNOTATOR_SCATTER.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH $PED_FILE \
# ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]} $CHROMOSOME $PHASE3_1KG_AUTOSOMES
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq);
# do
# CREATE_FAMILY_INFO_ARRAY
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_VARIANT_ANNOTATOR
# 		echo sleep 1s
# 	done
# done

# ##############################################################################################
# ##### GATHER UP THE PER FAMILY PER CHROMOSOME ANNOTATED VCF FILES INTO A SINGLE VCF FILE #####
# ##############################################################################################

# BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01_VARIANT_ANNOTATOR_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# CALL_VARIANT_ANNOTATOR_GATHER ()
# {
# echo \
# qsub \
# -N P.01-A.01_VARIANT_ANNOTATOR_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.ADD_MORE_ANNOTATION_GATHER.log \
#  $SCRIPT_DIR/P.01-A.01_VARIANT_ANNOTATOR_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }


# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
# 	BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_VARIANT_ANNOTATOR_GATHER
# 	echo sleep 1s
#  done

# ############################################################################################################
# ##### DO PER CHROMOSOME VARIANT TO TABLE FOR COHORT ########################################################
# ############################################################################################################

# CREATE_FAMILY_ONLY_ARRAY ()
# {
# FAMILY_ONLY_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $1,$20,$12,$18,$17}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_VARIANT_TO_TABLE_COHORT_ALL_SITES ()
# {
# echo \
# qsub \
# -N P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -hold_jid P.01_VARIANT_ANNOTATOR_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.VARIANT_TO_TABLE_COHORT_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
# do
# CREATE_FAMILY_ONLY_ARRAY
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_VARIANT_TO_TABLE_COHORT_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done

# ################################################################################################################
# ##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR COHORT ########################################################
# ################################################################################################################

# BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_COHORT_GATHER ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.02_VARIANT_TO_TABLE_COHORT_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# CALL_VARIANT_TO_TABLE_COHORT_GATHER ()
# {
# echo \
# qsub \
# -N T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
# 	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_COHORT_GATHER
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_VARIANT_TO_TABLE_COHORT_GATHER
# 	echo sleep 1s
#  done

# ##############################################################################################################
# ## BGZIP INITIAL JOINT CALLED VCF TABLE ######################################################################
# ##############################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$2"_"$1,\
# "-hold_jid","T.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# ##############################################################################################################
# ## TABIX INDEX INITIAL JOINT CALLED VCF TABLE ################################################################
# ##############################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","T.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES_"$2"_"$1,\
# "-hold_jid","T.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'


# #################################################################################
# ########### RUNNING FILTER TO FAMILY ALL SITES BY CHROMOSOME ####################
# #################################################################################

# CALL_FILTER_TO_FAMILY_ALL_SITES ()
# {
# echo \
# qsub \
# -N P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -hold_jid P.01_VARIANT_ANNOTATOR_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.FILTER_TO_FAMILY_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
# do
# CREATE_FAMILY_ONLY_ARRAY
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_FILTER_TO_FAMILY_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done
	
# #####################################################################################################
# ##### GATHER UP THE PER FAMILY PER CHROMOSOME FILTER TO FAMILY VCF FILES INTO A SINGLE VCF FILE #####
# #####################################################################################################

# BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# CALL_FILTER_TO_FAMILY_VCF_GATHER ()
# {
# echo \
# qsub \
# -N T.03-1_FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.FILTER_TO_FAMILY_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.03-1_FILTER_TO_FAMILY_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
# 	BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_FILTER_TO_FAMILY_VCF_GATHER
# 	echo sleep 1s
#  done

# ############################################################################################################
# ##### DO PER CHROMOSOME VARIANT TO TABLE FOR FAMILY ########################################################
# ############################################################################################################

# CALL_VARIANT_TO_TABLE_FAMILY_ALL_SITES ()
# {
# echo \
# qsub \
# -N T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -hold_jid P.01-A.03_FILTER_TO_FAMILY_ALL_SITES_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.VARIANT_TO_TABLE_FAMILY_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} $CHROMOSOME
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
# do
# CREATE_FAMILY_ONLY_ARRAY
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_VARIANT_TO_TABLE_FAMILY_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done

# ################################################################################################################
# ##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR FAMILY ########################################################
# ################################################################################################################

# BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_FAMILY_GATHER ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"T.03-2_VARIANT_TO_TABLE_FAMILY_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }


# CALL_VARIANT_TO_TABLE_FAMILY_GATHER ()
# {
# echo \
# qsub \
# -N T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.VARIANT_TO_TABLE_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
# 	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_FAMILY_GATHER
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_VARIANT_TO_TABLE_FAMILY_GATHER
# 	echo sleep 1s
#  done

# ##############################################################################################################
# ## BGZIP FAMILY ONLY VCF TABLE ###############################################################################
# ##############################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES_"$2"_"$1,\
# "-hold_jid","T.03-2-A.01_VARIANT_TO_TABLE_FAMILY_ALL_SITES_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# ##############################################################################################################
# ## TABIX INDEX FAMILY ONLY VCF TABLE #########################################################################
# ##############################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","T.03-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES_"$2"_"$1,\
# "-hold_jid","T.03-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_FAMILY_ALL_SITES_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.03-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_FAMILY_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# #################################################################################
# ########### RUNNING FILTER TO SAMPLE ALL SITES BY CHROMOSOME ####################
# #################################################################################

# # CREATE_SAMPLE_INFO_ARRAY_2 ()
# # {
# # SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {split($8,smtag,"[@]"); print $1,$8,$20,$12,$18,smtag[1]"_"smtag[2]}' \
# # ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# # for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
# # do
# # CREATE_SAMPLE_INFO_ARRAY_2
# # 	for CHROMOSOME in {{1..22},{X,Y}}
# # 		do
# # 		CALL_FILTER_TO_SAMPLE_ALL_SITES
# # 		echo sleep 1s
# # 		done
# # 	done

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_FILTER_TO_SAMPLE_ALL_SITES ()
# {
# echo \
# qsub \
# -N P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -hold_jid P.01_VARIANT_ANNOTATOR_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_2
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_FILTER_TO_SAMPLE_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done

# #####################################################################################################
# ##### GATHER UP THE PER SAMPLE PER CHROMOSOME FILTER TO SAMPLE VCF FILES INTO A SINGLE VCF FILE #####
# #####################################################################################################

# BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# # CALL_FILTER_TO_SAMPLE_VCF_GATHER ()
# # {
# # echo \
# # qsub \
# # -N T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]} \
# #  ${HOLD_ID_PATH} \
# #  -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[2]}/${SAMPLE_INFO_ARRAY_2[1]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_GATHER.log \
# #  $SCRIPT_DIR/T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER.sh \
# #  $JAVA_1_8 $GATK_DIR $CORE_PATH \
# #  ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[3]}
# # }

# # SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# # print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}'

# CALL_FILTER_TO_SAMPLE_VCF_GATHER ()
# {
# echo \
# qsub \
# -N T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.06-1_FILTER_TO_SAMPLE_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
# }

# # for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
# #  do
# #  	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF
# # 	CREATE_SAMPLE_INFO_ARRAY_2
# # 	CALL_FILTER_TO_SAMPLE_VCF_GATHER
# # 	echo sleep 1s
# #  done

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
#  do
#  	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF
# 	CREATE_SAMPLE_INFO_ARRAY_2
# 	CALL_FILTER_TO_SAMPLE_VCF_GATHER
# 	echo sleep 1s
#  done

# ############################################################################################################
# ##### DO PER CHROMOSOME VARIANT TO TABLE FOR SAMPLE ########################################################
# ############################################################################################################

# # CREATE_SAMPLE_INFO_ARRAY_2 ()
# # {
# # SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# # print $1,$20,$8,$12,$16,smtag[1]"_"smtag[2]}' \
# # ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# CALL_VARIANT_TO_TABLE_SAMPLE_ALL_SITES ()
# {
# echo \
# qsub \
# -N T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -hold_jid P.01-A.04_FILTER_TO_SAMPLE_ALL_SITES_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.VARIANT_TO_TABLE_SAMPLE_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_2
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_VARIANT_TO_TABLE_SAMPLE_ALL_SITES
# 		echo sleep 1s
# 		done
# 	done

# ################################################################################################################
# ##### GATHER PER CHROMOSOME VARIANT TO TABLE FOR SAMPLE ########################################################
# ################################################################################################################

# BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_SAMPLE_GATHER ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"T.06-2_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# CALL_VARIANT_TO_TABLE_SAMPLE_GATHER ()
# {
# echo \
# qsub \
# -N T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER.log \
#  $SCRIPT_DIR/T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
#  do
# 	BUILD_HOLD_ID_PATH_VARIANT_TO_TABLE_SAMPLE_GATHER
# 	CREATE_SAMPLE_INFO_ARRAY_2
# 	CALL_VARIANT_TO_TABLE_SAMPLE_GATHER
# 	echo sleep 1s
#  done

# #################################################################################################################
# ## ## BGZIP SAMPLE ONLY VCF TABLE ###############################################################################
# #################################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","T.06-2-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_GATHER_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# #################################################################################################################
# ## ## TABIX INDEX SAMPLE ONLY VCF TABLE #########################################################################
# #################################################################################################################

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","T.06-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","T.06-2-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.log",\
# "'$SCRIPT_DIR'""/T.06-2-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.sh",\
# "'$TABIX_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ###########################################################################################
# ########### RUNNING FILTER TO SAMPLE ALL SITES BY CHROMOSOME ON TARGET ####################
# ###########################################################################################

# # CREATE_SAMPLE_INFO_ARRAY_2 ()
# # {
# # SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$8,$20,$12,$18}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# CREATE_SAMPLE_INFO_ARRAY_2 ()
# {
# SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# print $1,$20,$8,$12,$18,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_FILTER_TO_SAMPLE_ALL_SITES_ON_TARGET ()
# {
# echo \
# qsub \
# -N P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -hold_jid P.01_VARIANT_ANNOTATOR_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]}_$CHROMOSOME \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_TARGET_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]} $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_2
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_FILTER_TO_SAMPLE_ALL_SITES_ON_TARGET
# 		echo sleep 1s
# 		done
# 	done

# ###############################################################################################################
# ##### GATHER UP THE PER SAMPLE PER CHROMOSOME FILTER TO SAMPLE VCF FILES ON TARGET INTO A SINGLE VCF FILE #####
# ###############################################################################################################

# BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF_TARGET ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.05_FILTER_TO_SAMPLE_ALL_SITES_TARGET_"$SAMPLE"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# # CREATE_SAMPLE_INFO_ARRAY_2 ()
# # {
# # SAMPLE_INFO_ARRAY_2=(`awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); if (smtag[1]"_"smtag[2]=="'$SAMPLE'") \
# # print $1,$20,$8,$12,$18,smtag[1]"_"smtag[2]}' \
# # ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# CALL_FILTER_TO_SAMPLE_VCF_TARGET_GATHER ()
# {
# echo \
# qsub \
# -N T.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER_${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE}_${SAMPLE_INFO_ARRAY_2[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${SAMPLE_INFO_ARRAY_2[0]}/${SAMPLE_INFO_ARRAY_2[1]}/${SAMPLE_INFO_ARRAY_2[2]}/LOGS/${SAMPLE_INFO_ARRAY_2[1]}_${SAMPLE_INFO_ARRAY_2[2]}_${SAMPLE_INFO_ARRAY_2[0]}.FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER.log \
#  $SCRIPT_DIR/T.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${SAMPLE_INFO_ARRAY_2[0]} ${SAMPLE_INFO_ARRAY_2[1]} ${SAMPLE_INFO_ARRAY_2[2]} ${SAMPLE_INFO_ARRAY_2[3]}
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {if ($8~"@") {split($8,smtag,"[@]"); print smtag[1]"_"smtag[2]} else print $8"_"}' $SAMPLE_SHEET | sort | uniq );
#  do
#  	BUILD_HOLD_ID_PATH_FILTER_TO_SAMPLE_VCF_TARGET
# 	CREATE_SAMPLE_INFO_ARRAY_2
# 	CALL_FILTER_TO_SAMPLE_VCF_TARGET_GATHER
# 	echo sleep 1s
#  done

# ###########################################################################################
# ########### RUNNING FILTER TO FAMILY ALL SITES BY CHROMOSOME ON TARGET ####################
# ###########################################################################################

# # CREATE_FAMILY_ONLY_ARRAY ()
# # {
# # FAMILY_ONLY_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $20=="'$FAMILY'" {print $1,$20,$12,$18,$17}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# # }

# CALL_FILTER_TO_FAMILY_ON_TARGET_VARIANT ()
# {
# echo \
# qsub \
# -N P.01-A.06_FILTER_TO_FAMILY_TARGET_VARIANT_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -hold_jid P.01_VARIANT_ANNOTATOR_${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}_$CHROMOSOME \
# -o $CORE_PATH/${FAMILY_ONLY_ARRAY[0]}/${FAMILY_ONLY_ARRAY[1]}/LOGS/${FAMILY_ONLY_ARRAY[1]}_${FAMILY_ONLY_ARRAY[0]}.FILTER_TO_FAMILY_ALL_SITES_$CHROMOSOME.log \
# $SCRIPT_DIR/P.01-A.06_FILTER_TO_FAMILY_ON_TARGET_VARIANT_ONLY_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH \
# ${FAMILY_ONLY_ARRAY[0]} ${FAMILY_ONLY_ARRAY[1]} ${FAMILY_ONLY_ARRAY[2]} ${FAMILY_ONLY_ARRAY[4]} $CHROMOSOME
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"} {print $1}' $PED_FILE | sort | uniq );
# do
# CREATE_FAMILY_ONLY_ARRAY
# 	for CHROMOSOME in {{1..22},{X,Y}}
# 		do
# 		CALL_FILTER_TO_FAMILY_ON_TARGET_VARIANT
# 		echo sleep 1s
# 		done
# 	done
	
# ###############################################################################################################
# ##### GATHER UP THE PER FAMILY PER CHROMOSOME ON TARGET FILTER TO FAMILY VCF FILES INTO A SINGLE VCF FILE #####
# ###############################################################################################################

# BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF_TARGET_VARIANT ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"P.01-A.06_FILTER_TO_FAMILY_TARGET_VARIANT_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
#  	done
#  done
# }

# CALL_FILTER_TO_FAMILY_VCF_GATHER_TARGET_VARIANT ()
# {
# echo \
# qsub \
# -N T.09-1_FILTER_TO_FAMILY_ON_TARGET_VARIANT_GATHER_${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]} \
#  ${HOLD_ID_PATH} \
#  -o $CORE_PATH/${FAMILY_INFO_ARRAY[0]}/${FAMILY_INFO_ARRAY[2]}/LOGS/${FAMILY_INFO_ARRAY[2]}_${FAMILY_INFO_ARRAY[0]}.FILTER_TO_FAMILY_ON_TARGET_VARIANT_GATHER.log \
#  $SCRIPT_DIR/T.09-1_FILTER_TO_FAMILY_ON_TARGET_VARIANT_ONLY_GATHER.sh \
#  $JAVA_1_8 $GATK_DIR $CORE_PATH \
#  ${FAMILY_INFO_ARRAY[0]} ${FAMILY_INFO_ARRAY[2]} ${FAMILY_INFO_ARRAY[3]}
# }

# for FAMILY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt | sort | uniq)
#  do
# 	BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF_TARGET_VARIANT
# 	CREATE_FAMILY_INFO_ARRAY
# 	CALL_FILTER_TO_FAMILY_VCF_GATHER_TARGET_VARIANT
# 	echo sleep 1s
#  done

# ###############################
# ##### DOING VCF BREAKOUTS #####
# ###############################

# ### SUBSETTING FROM COHORT (FAMILY PLUS CONTROL SET) VCF ###

# # FILTER TO JUST VARIANT SITES
# # I think Molly might like this output, but if not, then don't have to generate it.

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","S.01_FILTER_COHORT_VARIANT_ONLY_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY.log",\
# "'$SCRIPT_DIR'""/S.01_FILTER_COHORT_VARIANT_ONLY.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # FILTER TO JUST PASSING VARIANT SITES
# # I think statgen is using this for some of their programs
# # If not needed then don't generate

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{print "qsub","-N","S.02_FILTER_COHORT_VARIANT_ONLY_PASS_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_VARIANT_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.02_FILTER_COHORT_VARIANT_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # FILTER TO JUST PASSING BIALLELIC SNV SITES
# # TEMPORARY FILE USED FOR PCA AND RELATEDNESS

# awk 'BEGIN {OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{print "qsub","-N","S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC.log",\
# "'$SCRIPT_DIR'""/S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # RUN HUAS WORKFLOW FOR PCA AND RELATEDNESS

# awk 'BEGIN {OFS="\t"} {print $1,$20,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{print "qsub","-N","S.03-A.01_PCA_RELATEDNESS_"$2"_"$1,\
# "-hold_jid","S.03_FILTER_COHORT_SNV_ONLY_PASS_BIALLELIC_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/LOGS/"$2"_"$1".PCA_RELATEDNESS.log",\
# "'$SCRIPT_DIR'""/S.03-A.01_PCA_RELATEDNESS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VCFTOOLS_DIR'","'$PLINK2_DIR'","'$KING_DIR'",$1,$2,$3,"'$PED_FILE'","'$CONTROL_PED_FILE'""\n""sleep 1s"}'

# #################################
# ### SUBSETTING TO SAMPLE VCFS ###
# #################################

# ## SUBSET TO SAMPLE VARIANTS ONLY ON BAIT

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.07_FILTER_TO_SAMPLE_VARIANTS_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".FILTER_TO_VARIANTS.log",\
# "'$SCRIPT_DIR'""/S.07_FILTER_TO_SAMPLE_VARIANTS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

# ## SUBSET TO SAMPLE PASSING SNVS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_SNV_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.09_FILTER_TO_SAMPLE_SNV_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO SAMPLE PASSING INDELS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.10_FILTER_TO_INDEL_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_INDEL_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.10_FILTER_TO_SAMPLE_INDEL_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO SAMPLE PASSING MIXED

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.11_FILTER_TO_MIXED_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_MIXED_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.11_FILTER_TO_SAMPLE_MIXED_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO TARGET SNV ONLY PASS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_SNV_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO TARGET INDEL ONLY PASS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_INDEL_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO TARGET MIXED ONLY PASS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TARGET_MIXED_ONLY_PASS.log",\
# "'$SCRIPT_DIR'""/S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# ## SUBSET TO SAMPLE VARIANTS ONLY ON TARGET

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","P.01-A.01_VARIANT_ANNOTATOR_GATHER_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".FILTER_TO_VARIANTS_TARGET.log",\
# "'$SCRIPT_DIR'""/S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'


# ####################
# ### TITV SECTION ###
# ####################

# # BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF.log",\
# "'$SCRIPT_DIR'""/S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# # BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND OVERLAP WITH DBSNP 129

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_KNOWN.log",\
# "'$SCRIPT_DIR'""/S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,"'$DBSNP_129'""\n""sleep 1s"}'

# # BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND DO NOT OVERLAP WITH DBSNP 129

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".FILTER_TO_TITV_VCF_NOVEL.log",\
# "'$SCRIPT_DIR'""/S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,"'$DBSNP_129'""\n""sleep 1s"}'

# ### RUN TITV FOR THE PASSING SNVS THAT FALL IN UCSC CODING REGIONS THAT TOUCH EITHER THE BED OR TARGET FILE

# ## ALL SNVS TITV

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.01-A.01_TITV_ALL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_ALL.log",\
# "'$SCRIPT_DIR'""/S.09-A.01-A.01_TITV_ALL.sh",\
# "'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ## ALL KNOWN SNVS TITV

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.02-A.01_TITV_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_KNOWN.log",\
# "'$SCRIPT_DIR'""/S.09-A.02-A.01_TITV_KNOWN.sh",\
# "'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ## ALL NOVEL SNVS TITV

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.09-A.03-A.01_TITV_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_TITV_NOVEL.log",\
# "'$SCRIPT_DIR'""/S.09-A.03-A.01_TITV_NOVEL.sh",\
# "'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# ###################
# ##### ANNOVAR #####
# ###################

# ## RUN ANNOVAR

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 -k 3,3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16-A.01_RUN_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.16_FILTER_TO_SAMPLE_VARIANTS_TARGET_"smtag[1]"_"smtag[2]"_"$1,\
# "-pe slots 5",\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".RUN_ANNOVAR.log",\
# "'$SCRIPT_DIR'""/S.16-A.01_RUN_ANNOVAR.sh",\
# "'$JAVA_1_6'","'$CIDRSEQSUITE_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

# ## REFORMAT ANNOVAR

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 -k 3,3 \
# | uniq \
# | awk '{split($3,smtag,"[@]"); print "qsub","-N","S.16-A.01-A.01_REFORMAT_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-hold_jid","S.16-A.01_RUN_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/"$2"/"$3"/LOGS/"$3"_"$2"_"$1".REFORMAT_ANNOVAR.log",\
# "'$SCRIPT_DIR'""/S.16-A.01-A.01_REFORMAT_ANNOVAR.sh",\
# "'$ANNOVAR_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

# ######### FINISH UP #################

# ### QC REPORT PREP ###

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$21,$22,$23,$24}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 -k 3 \
# | uniq \
# | awk 'BEGIN {FS="\t"}
# {split($3,smtag,"[@]"); print "qsub","-N","X.01-QC_REPORT_PREP_"$1"_"smtag[1]"_"smtag[2],\
# "-hold_jid","S.16-A.01-A.01_REFORMAT_ANNOVAR_"smtag[1]"_"smtag[2]"_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".QC_REPORT_PREP.log",\
# "'$SCRIPT_DIR'""/X.01-QC_REPORT_PREP.sh",\
# "'$SAMTOOLS_DIR'","'$CORE_PATH'","'$DATAMASH_DIR'",$1,$2,$3,$4,$5,$6,$7"\n""sleep 1s"}'

# ### END PROJECT TASKS ###

# awk 'BEGIN {FS="\t"; OFS="\t"} {split($8,smtag,"[@]"); print $1,smtag[1]"_"smtag[2]}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | $DATAMASH_DIR/datamash -s -g 1 collapse 2 \
# | awk 'BEGIN {FS="\t"}
# gsub (/,/,",X.01-QC_REPORT_PREP_"$1"_",$2) \
# {print "qsub","-N","X.01-X.01-END_PROJECT_TASKS_"$1,\
# "-hold_jid","X.01-QC_REPORT_PREP_"$1"_"$2,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$1".END_PROJECT_TASKS.log",\
# "'$SCRIPT_DIR'""/X.01-X.01-END_PROJECT_TASKS.sh",\
# "'$CORE_PATH'","'$DATAMASH_DIR'",$1"\n""sleep 1s"}'
