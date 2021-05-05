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

	THREADS=$6 # optional. if no 6th argument present then default is 6.
		# if you wangt to set this then you need to set 3rd,4th and 5th argument as well (even to default)

			if [[ ! $THREADS ]]
				then
				THREADS="6"
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

	# grab email addy

		SEND_TO=`cat $SCRIPT_DIR/../email_lists.txt`

	# grab submitter's name

		PERSON_NAME=`getent passwd | awk 'BEGIN {FS=":"} $1=="'$SUBMITTER_ID'" {print $5}'`

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

	ALIGNMENT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/ddl_ce_control_align-0.0.4.simg"
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

	GATK_3_7_0_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/gatk3-3.7-0.simg"
	# singularity pull docker://broadinstitute/gatk3:3.7-0
	# used for generating the depth of coverage reports.
		# comes with R 3.1.1 with appropriate packages needed to create gatk pdf output
		# also comes with some version of java 1.8
		# jar file is /usr/GenomeAnalysisTK.jar

	MITO_MUTECT2_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/mito_mutect2-4.1.3.0.0.simg"
		# uses broadinstitute/gatk:4.1.3.0 as the base image (as /gatk/gatk.jar)
			# added
				# bcftools-1.10.2
				# haplogrep-2.1.20.jar (as /jars/haplogrep-2.1.20.jar)
				# annovar

	MITO_EKLIPSE_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/mito_eklipse-master-c25931b.0.simg"
		# https://github.com/dooguypapua/eKLIPse AND all of its dependencies

	MT_COVERAGE_R_SCRIPT="$SCRIPT_DIR/mito_coverage_graph.r"

	CNV_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/cnv_exomedepth-dev.simg"

	EXOME_DEPTH_R_SCRIPT="$SCRIPT_DIR/runExomeDepth.r"

	FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT="$SCRIPT_DIR/FORMAT_AND_ZOOM_ANNOTSV.r"

	PCA_RELATEDNESS_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/pca-relatedness-0.0.1.simg"

	VT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/vt-0.5772.ca352e2c.0.simg"

	ANNOVAR_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/annovarwrangler-20210126.simg"

	# PIPELINE PROGRAMS TO BE IMPLEMENTED
		# JAVA_1_6="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jre1.6.0_25/bin"
		# CIDRSEQSUITE_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/CIDRSeqSuiteSoftware_Version_4_0/"
		# ANNOVAR_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/ANNOVAR/2013_09_11"

	# ANNOVAR PARAMETERS AND INPUTS

		ANNOVAR_DATABASE_FILE="$SCRIPT_DIR/../resources/CFTR.final.csv"
		ANNOVAR_REF_BUILD="hg19"

		ANNOVAR_INFO_FIELD_KEYS="VariantType," \
			ANNOVAR_INFO_FIELD_KEYS=$ANNOVAR_INFO_FIELD_KEYS"DP" \

		ANNOVAR_HEADER_MAPPINGS="af=gnomad211_exome_AF," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_popmax=gnomad211_exome_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_male=gnomad211_exome_AF_male," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_female=gnomad211_exome_AF_female," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_raw=gnomad211_exome_AF_raw," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_afr=gnomad211_exome_AF_afr," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_sas=gnomad211_exome_AF_sas," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_amr=gnomad211_exome_AF_amr," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_eas=gnomad211_exome_AF_eas," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_nfe=gnomad211_exome_AF_nfe," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_fin=gnomad211_exome_AF_fin," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_asj=gnomad211_exome_AF_asj," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"af_oth=gnomad211_exome_AF_oth," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_topmed_af_popmax=gnomad211_exome_non_topmed_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_neuro_af_popmax=gnomad211_exome_non_neuro_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_cancer_af_popmax=gnomad211_exome_non_cancer_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"controls_af_popmax=gnomad211_exome_controls_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF=gnomad211_genome_AF," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_popmax=gnomad211_genome_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_male=gnomad211_genome_AF_male," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_female=gnomad211_genome_AF_female," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_raw=gnomad211_genome_AF_raw," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_afr=gnomad211_genome_AF_afr," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_sas=gnomad211_genome_AF_sas," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_amr=gnomad211_genome_AF_amr," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_eas=gnomad211_genome_AF_eas," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_nfe=gnomad211_genome_AF_nfe," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_fin=gnomad211_genome_AF_fin," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_asj=gnomad211_genome_AF_asj," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"AF_oth=gnomad211_genome_AF_oth," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_topmed_AF_popmax=gnomad211_genome_non_topmed_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_neuro_AF_popmax=gnomad211_genome_non_neuro_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"non_cancer_AF_popmax=gnomad211_genome_non_cancer_AF_popmax," \
			ANNOVAR_HEADER_MAPPINGS=$ANNOVAR_HEADER_MAPPINGS"controls_AF_popmax=gnomad211_genome_controls_AF_popmax"

			ANNOVAR_VCF_COLUMNS="CHROM,"
				ANNOVAR_VCF_COLUMNS=$ANNOVAR_VCF_COLUMNS"POS,"
				ANNOVAR_VCF_COLUMNS=$ANNOVAR_VCF_COLUMNS"REF,"
				ANNOVAR_VCF_COLUMNS=$ANNOVAR_VCF_COLUMNS"ALT"

##################
# PIPELINE FILES #
##################

	# Core Pipeline

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

		# where the control data set resides.

		CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Data/TWIST_CONTROL_SET1.200601_PIPELINE_2_0_0"
		CONTROL_PED_FILE="$CONTROL_REPO/TWIST_CONTROL_SET1.200601.ped"

		# SFAFASFA

		CONTROL_DATA_SET_FILE="CGC_CONTROL_SET_3_7.g.vcf.gz"

	# Mitochondrial Pipline

		MT_PICARD_INTERVAL_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/MT.interval_list"
		MT_MASK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/hg37_MT_blacklist_sites.hg37.MT.bed"
		GNOMAD_MT="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/GRCh37_MT_gnomAD.vcf.gz"
		ANNOVAR_MT_DB_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/annovar_db/"
		MT_GENBANK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/NC_012920.1.gb"

	# CNV calling workflow

		## REF_PANEL_COUNTS USED IN EXOME DEPTH IS SEX SPECIFIC.
		## DETERMINED WHEN PARSING GENDER FROM PED FILE DURING CREATE_SAMPLE_ARRAY
		## THE THREE RDA FILES BELOW GET REASSIGNED TO $REF_PANEL_COUNTS depending on what the gender is.

			# read count from female reference panel, won't need to change unless changes in bed file or reference samples

				REF_PANEL_FEMALE_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountFemaleUniqBed22.rda"

			# read count from male reference panel, won't need to change unless changes in bed file or reference samples

				REF_PANEL_MALE_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountMaleUniqBed26.rda"

			# if subject sex is not specified as 'm' or 'f', it will use count of all sample

				REF_PANEL_ALL_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountAllUniqBed48.rda"

#################################
##### MAKE A DIRECTORY TREE #####
#################################

	# make a directory in user home directory

		mkdir -p ~/CGC_PIPELINE_TEMP

	# create variables using the base name for the sample sheet and ped file

		MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
		PED_PREFIX=`basename $PED_FILE .ped`

	# fix any commonly seen formatting issues in the sample sheet

		FORMAT_MANIFEST ()
		{
			awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'NR>1' \
			| sed 's/,/\t/g' \
			| sort -k 8,8 \
			>| ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt
		}

	# merge the sample sheet with the ped file

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

	# create an array from values of the merged sample sheet and ped file
	## set $REF_PANEL_COUNTS based on gender

		CREATE_SAMPLE_ARRAY ()
		{
			SAMPLE_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" \
				{split($19,INDEL,";"); \
				print $1,$8,$9,$10,$11,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$20,$21,$22,$23,$24}' \
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

			# 10  Description=Sequencer model and/or setting (setting e.g. "Rapid-Run")
			## Models: “HiSeq-X”,“HiSeq-4000”,“HiSeq-2500”,“HiSeq-2000”,“NextSeq-500”,“MiSeq”

				SEQUENCER_MODEL=${SAMPLE_ARRAY[3]}

			# 11  Seq_Exp_ID=Zoom Gene List for filtering ExomeDepth output by gene symbol

				ZOOM_LIST=${SAMPLE_ARRAY[4]}

					# if the zoom list file exists than the output file prefix is the input file prefix before .GeneList

						if [ -f $ZOOM_LIST ]
							then ZOOM_NAME=$(basename $ZOOM_LIST | sed 's/.GeneList.[0-9]*.csv//g')
							else ZOOM_NAME="NA"
						fi

			# 12  Genome_Ref=the reference genome used in the analysis pipeline

				REF_GENOME=${SAMPLE_ARRAY[5]}
					REF_DICT=$(echo $REF_GENOME | sed 's/fasta$/dict/g; s/fa$/dict/g')

				#####################################
				# 13  Operator: SKIP ################
				# 14  Extra_VCF_Filter_Params: SKIP #
				#####################################

			# 15  TS_TV_BED_File=where ucsc coding exons overlap with bait and target bed files

				TITV_BED=${SAMPLE_ARRAY[6]}

			# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
			# Used for limited where to run base quality score recalibration on where to create gvcf files.

				BAIT_BED=${SAMPLE_ARRAY[7]}

			# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.

				TARGET_BED=${SAMPLE_ARRAY[8]}

			# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.

				DBSNP=${SAMPLE_ARRAY[9]}

			# 19  KNOWN_INDEL_FILES=used for BQSR masking, sensitivity in local realignment.

				KNOWN_INDEL_1=${SAMPLE_ARRAY[10]}
				KNOWN_INDEL_2=${SAMPLE_ARRAY[11]}

			# 20 family that sample belongs to

				FAMILY=${SAMPLE_ARRAY[12]}

			# 21 MOM

				FATHER=${SAMPLE_ARRAY[13]}

			# 22 DAD

				MOTHER=${SAMPLE_ARRAY[14]}

			# 23 GENDER

				GENDER=${SAMPLE_ARRAY[15]}

				# set $REF_PANEL_COUNTS USED IN EXOMEDEPTH TO THE SEX SPECIFIC ONE

					if [[ $GENDER = "1" ]];
						then REF_PANEL_COUNTS=${REF_PANEL_MALE_READ_COUNT_RDA}
					elif [[ $GENDER = "2" ]];
						then REF_PANEL_COUNTS=${REF_PANEL_FEMALE_READ_COUNT_RDA}
					else
						REF_PANEL_COUNTS=${REF_PANEL_ALL_READ_COUNT_RDA}
					fi

			# 24 PHENOTYPE

				PHENOTYPE=${SAMPLE_ARRAY[16]}
		}

	# PROJECT DIRECTORY TREE CREATOR

		MAKE_PROJ_DIR_TREE ()
		{
			mkdir -p $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CNV_OUTPUT \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/CRAM \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/HC_CRAM \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/GVCF \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_AUTO,RG_HEADER,QUALITY_YIELD,ERROR_SUMMARY,VCF_METRICS} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/CONCORDANCE \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/DEPTH_OF_COVERAGE/{TARGET_PADDED,CODING_PADDED} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/INSERT_SIZE/{METRICS,PDF} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/REPORTS/ANEUPLOIDY_CHECK \
			$CORE_PATH/$PROJECT/$FAMILY/{LOGS,VCF,RELATEDNESS,PCA} \
			$CORE_PATH/$PROJECT/$FAMILY/VCF/{RAW,VQSR} \
			$CORE_PATH/$PROJECT/TEMP/${SM_TAG}_ANNOVAR_TARGET \
			$CORE_PATH/$PROJECT/TEMP/{VCF_PREP,PLINK,KING} \
			$CORE_PATH/$PROJECT/{TEMP,FASTQ,REPORTS,LOGS,COMMAND_LINES} \
			$CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/MT_OUTPUT/{COLLECTHSMETRICS_MT,MUTECT2_MT,HAPLOTYPES,ANNOVAR_MT,EKLIPSE} \
			$CORE_PATH/$PROJECT/TEMP/${SM_TAG}_ANNOVAR_MT
		}

	# combine above functions into one...this is probably not necessary...

		SETUP_PROJECT ()
		{
			FORMAT_MANIFEST
			MERGE_PED_MANIFEST
			CREATE_SAMPLE_ARRAY
			MAKE_PROJ_DIR_TREE
			echo Project started at `date` >| $CORE_PATH/$PROJECT/REPORTS/PROJECT_START_END_TIMESTAMP.txt
		}

############################################
# run steps for pipeline and project setup #
############################################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		SETUP_PROJECT
done

#######################################################################################
##### CRAM FILE GENERATION ############################################################
# NOTE: THE CRAM FILE IS THE END PRODUCT BUT THE BAM FILE IS USED FOR OTHER PROCESSES #
# SOME PROGRAMS CAN'T TAKE IN CRAM AS AN INPUT ########################################
#######################################################################################

	########################################################################################
	# create an array at the platform level so that bwa mem can add metadata to the header #
	########################################################################################

		CREATE_PLATFORM_UNIT_ARRAY ()
		{
			PLATFORM_UNIT_ARRAY=(`awk 1 ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS="\t"} $8$2$3$4=="'$PLATFORM_UNIT'" {split($19,INDEL,";"); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$20,$21,$22,$23,$24}' \
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

				# 10  Description=Sequencer model and/or setting (setting e.g. "Rapid-Run")
				## Models: “HiSeq-X”,“HiSeq-4000”,“HiSeq-2500”,“HiSeq-2000”,“NextSeq-500”,“MiSeq”

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

					MOM=${PLATFORM_UNIT_ARRAY[18]}

				# 22 DAD

					DAD=${PLATFORM_UNIT_ARRAY[19]}

				# 23 GENDER

					GENDER=${PLATFORM_UNIT_ARRAY[20]}

				# 24 PHENOTYPE

					PHENOTYPE=${PLATFORM_UNIT_ARRAY[21]}
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
				$NOVASEQ_REPO \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	for PLATFORM_UNIT in $(awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $8$2$3$4}' \
			| sort \
			| uniq );
		do
			CREATE_PLATFORM_UNIT_ARRAY
			mkdir -p $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/
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
				"'$THREADS'",\
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

	#######################################
	# run bqsr on the using bait bed file #
	#######################################

		PERFORM_BQSR ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-PERFORM_BQSR.log" \
			-hold_jid B.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT,C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/D.01_PERFORM_BQSR.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$KNOWN_INDEL_1 \
				$KNOWN_INDEL_2 \
				$DBSNP \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################
	# use a 4 bin q score scheme #
	# remove indel Q scores ######
	# retain original Q score  ###
	##############################

		APPLY_BQSR ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-APPLY_BQSR.log" \
			-hold_jid D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/E.01_APPLY_BQSR.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#####################################################
	# create a lossless cram, although the bam is lossy #
	#####################################################

		BAM_TO_CRAM ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-BAM_TO_CRAM.log" \
			-hold_jid E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/F.01_BAM_TO_CRAM.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

######################################
# run steps for cram file generation #
######################################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		CREATE_SAMPLE_ARRAY
		FIX_BED_FILES
		echo sleep 0.1s
		PERFORM_BQSR
		echo sleep 0.1s
		APPLY_BQSR
		echo sleep 0.1s
		BAM_TO_CRAM
		echo sleep 0.1s
done

####################################################################################
##### BAM/CRAM FILE RELATED METRICS ################################################
# NOTE: SOME PROGRAMS CAN ONLY BE RAN ON THE BAM FILE AND NOT ON THE CRAM FILE #####
# I WILL COMMENT ON WHICH IS WHICH #################################################
####################################################################################

	################################################################################
	# COLLECT MULTIPLE METRICS  ####################################################
	# again used bait bed file here instead of target b/c target could be anything #
	# ti/tv bed is unrelated to the capture really #################################
	# uses the CRAM file as the input ##############################################
	################################################################################

		COLLECT_MULTIPLE_METRICS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.01-COLLECT_MULTIPLE_METRICS"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-COLLECT_MULTIPLE_METRICS.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.01_COLLECT_MULTIPLE_METRICS.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$DBSNP \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#########################################
	# COLLECT HS METRICS  ###################
	# bait bed is the bait bed file #########
	# titv bed files is the target bed file #
	# uses the CRAM file as the input #######
	#########################################

		COLLECT_HS_METRICS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.02-COLLECT_HS_METRICS"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-COLLECT_HS_METRICS.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.02_COLLECT_HS_METRICS.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$BAIT_BED \
				$TITV_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR TARGET BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container ##################################################
	# input is the BAM file #################################################################################
	# Generally this with all RefSeq Select CDS exons + missing OMIM unless it becomes targeted, e.g a zoom #
	# uses the BAM file as the input ########################################################################
	#########################################################################################################

		DOC_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.03-DOC_TARGET"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-DOC_TARGET.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.03_DOC_TARGET_PADDED_BED.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$TARGET_BED \
				$PADDING_LENGTH \
				$GENE_LIST \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#################################################################################################
	# CREATE VCF FOR VERIFYBAMID METRICS ############################################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	#################################################################################################

		SELECT_VERIFYBAMID_VCF ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.04-SELECT_VERIFYBAMID_VCF"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-SELECT_VERIFYBAMID_VCF.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.04_SELECT_VERIFYBAMID_VCF.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$VERIFY_VCF \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	###############################
	# RUN VERIFYBAMID #############
	# THIS RUNS OFF OF A BAM FILE #
	###############################

		RUN_VERIFYBAMID ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.04-A.01-RUN_VERIFYBAMID"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-VERIFYBAMID.log" \
			-hold_jid H.04-SELECT_VERIFYBAMID_VCF"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.04-A.01_VERIFYBAMID.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR CODING BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container ##################################################
	# input is the BAM file ######################################################
	# This with all RefSeq Select CDS exons + missing OMIM, etc. #################
	# uses the BAM file as the input #############################################
	##############################################################################

		DOC_CODING ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-DOC_CODING.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05_DOC_CODING_PADDED.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$CODING_BED \
				$PADDING_LENGTH \
				$GENE_LIST \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#########################################################
	# DO AN ANEUPLOIDY CHECK ON TARGET BED FILE DOC OUTPUT  #
	#########################################################

		ANEUPLOIDY_CHECK ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.01_CHROM_DEPTH"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-ANEUPLOIDY_CHECK.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.01_CHROM_DEPTH.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	########################################################################################
	# FORMATTING PER BASE COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	########################################################################################

		ANNOTATE_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-ANNOTATE_PER_BASE.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02_ANNOTATE_PER_BASE.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##########################################################################
	# FILTER PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x #
	##########################################################################

		FILTER_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02-A.01_FILTER_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-FILTER_ANNOTATED_PER_BASE.log" \
			-hold_jid H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02-A.01_FILTER_ANNOTATED_PER_BASE.sh \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	######################################################
	# BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION #
	######################################################

		BGZIP_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02-A.02_BGZIP_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-BGZIP_ANNOTATED_PER_BASE.log" \
			-hold_jid H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02-A.02_BGZIP_ANNOTATED_PER_BASE.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	###################################################################################################
	# FORMATTING PER CODING INTERVAL COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	###################################################################################################

		ANNOTATE_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.03_ANNOTATE_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-ANNOTATE_PER_INTERVAL.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.03_ANNOTATE_PER_INTERVAL.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	##################################################################################################
	# FILTER ANNOTATED PER CODING INTERVAL COVERAGE TO INTERVALS WHERE LESS 100% OF BASES ARE AT 30X #
	##################################################################################################

		FILTER_ANNOTATED_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.03-A.01_FILTER_ANNOTATED_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-FILTER_ANNOTATED_PER_INTERVAL.log" \
			-hold_jid H.05-A.03_ANNOTATE_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.03-A.01_FILTER_ANNOTATED_PER_INTERVAL.sh \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	#################################################################################################
	# CREATE VCF PER CHROMOSOME AND RUN VERIFYBAMID ON THEM ######################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	# SCRIPT READS BAIT BED FILE, GRABS THE CHROMOSOMES AND RUNS A FOR LOOP FOR BOTH THINGS #########
	# USES BAM FILE AS THE INPUT ####################################################################
	#################################################################################################

		VERIFYBAMID_PER_AUTOSOME ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.06-SELECT_VERIFYBAMID_PER_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-SELECT_VERIFYBAMID_PER_AUTOSOME.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.06_VERIFYBAMID_PER_AUTO.sh \
				$ALIGNMENT_CONTAINER \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$VERIFY_VCF \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#################################################################################################
	# CREATE VCF PER CHROMOSOME AND RUN VERIFYBAMID ON THEM ######################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	# SCRIPT READS BAIT BED FILE, GRABS THE CHROMOSOMES AND RUNS A FOR LOOP FOR BOTH THINGS #########
	# USES BAM FILE AS THE INPUT ####################################################################
	#################################################################################################

		CAT_VERIFYBAMID_PER_AUTOSOME ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.06-A.01-CAT_VERIFYBAMID_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-CAT_VERIFYBAMID_AUTOSOME.log" \
			-hold_jid H.06-SELECT_VERIFYBAMID_PER_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.06-A.01_CAT_VERIFYBAMID_AUTO.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$BAIT_BED
		}

###############################################
# run steps for cram/bam file related metrics #
###############################################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $8}' \
			| sort \
			| uniq );
	do
		CREATE_SAMPLE_ARRAY
		COLLECT_MULTIPLE_METRICS
		echo sleep 0.1s
		COLLECT_HS_METRICS
		echo sleep 0.1s
		DOC_TARGET
		echo sleep 0.1s
		SELECT_VERIFYBAMID_VCF
		echo sleep 0.1s
		RUN_VERIFYBAMID
		echo sleep 0.1s
		DOC_CODING
		echo sleep 0.1s
		ANEUPLOIDY_CHECK
		echo sleep 0.1s
		ANNOTATE_PER_BASE_REPORT
		echo sleep 0.1s
		FILTER_ANNOTATED_PER_BASE_REPORT
		echo sleep 0.1s
		BGZIP_ANNOTATED_PER_BASE_REPORT
		echo sleep 0.1s
		ANNOTATE_PER_INTERVAL_REPORT
		echo sleep 0.1s
		FILTER_ANNOTATED_PER_INTERVAL_REPORT
		echo sleep 0.1s
		VERIFYBAMID_PER_AUTOSOME
		echo sleep 0.1s
		CAT_VERIFYBAMID_PER_AUTOSOME
		echo sleep 0.1s
done

#############################################
##### MITOCHONDRIAL WORKFLOW ################
# RUN MUTECT2 TO CALL SNVS AND SMALL INDELS #
# GENERATE COVERAGE STATS ###################
# RUN EKLIPSE FOR LARGER DELETIONS ##########
#############################################

	#########################################
	##### MUTECT2 IN MITO MODE WORKFLOW #####
	##### WORKS ON FULL BAM FILE ############
	#########################################

		#####################################################
		# run mutect2 in mitochondria mode on full bam file #
		# this runs MUCH slower on non-avx machines #########
		#####################################################

			MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-MUTECT2_MT.log" \
					-hold_jid E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$REF_GENOME \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		#######################################
		# apply filters to mutect2 vcf output #
		#######################################

			FILTER_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-FILTER_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-FILTER_MUTECT2_MT.log" \
					-hold_jid H.08-MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-FILTER_MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$SM_TAG \
					$REF_GENOME \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		###################################################
		# apply masks to mutect2 mito filtered vcf output #
		###################################################

			MASK_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-A.01-MASK_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-MASK_MUTECT2_MT.log" \
					-hold_jid H.08-A.01-FILTER_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-A.01-MASK_MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$SM_TAG \
					$MT_MASK \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		#############################################
		# run haplogrep2 on mutect2 mito vcf output #
		#############################################

			HAPLOGREP2_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-A.01-A01-HAPLOGREP2_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-HAPLOGREP2_MUTECT2_MT.log" \
					-hold_jid H.08-A.01-A.01-MASK_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-A.01-A.01-HAPLOGREP2_MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$REF_GENOME \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		#############################################################
		# add gnomad annotation to info field of mutect2 vcf output #
		#############################################################

			GNOMAD_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-A.01-A.02-GNOMAD_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-GNOMAD_MUTECT2_MT.log" \
					-hold_jid H.08-A.01-A.01-MASK_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-A.01-A.02-GNOMAD_MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$GNOMAD_MT \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		##########################################
		# run annovar on final mutect2 based vcf #
		##########################################

			RUN_ANNOVAR_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-A.01-A.02-A01-RUN_ANNOVAR_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-RUN_ANNOVAR_MUTECT2_MT.log" \
					-hold_jid H.08-A.01-A.01-A.02-GNOMAD_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-A.01-A.02-A.01-RUN_ANNOVAR_MUTECT2_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$ANNOVAR_MT_DB_DIR \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		##########################################
		# run annovar on final mutect2 based vcf #
		##########################################

			FIX_ANNOVAR_MUTECT2_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.08-A.01-A.01-A.02-A.01-A.01-FIX_ANNOVAR_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-FIX_ANNOVAR_MUTECT2_MT.log" \
					-hold_jid H.08-A.01-A.01-A.02-A01-RUN_ANNOVAR_MUTECT2_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.08-A.01-A.01-A.02-A.01-A.01-FIX_ANNOVAR_MUTECT2_MT.sh \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

	##############################################################
	##### RUN EKLIPSE TO DETECT LARGE DELETIONS IN MT GENOME #####
	##############################################################

		############################################
		# SUBSET BAM FILE TO CONTAIN ONLY MT READS #
		############################################

			SUBSET_BAM_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.09-MAKE_BAM_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-MAKE_BAM_MT.log" \
					-hold_jid E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.09-MAKE_MT_BAM.sh \
					$MITO_EKLIPSE_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$THREADS \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		###############
		# RUN EKLIPSE #
		###############

			RUN_EKLIPSE ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.09-A.01-RUN_EKLIPSE"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-RUN_EKLIPSE.log" \
					-hold_jid H.09-MAKE_BAM_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.09-A.01-RUN_EKLIPSE.sh \
					$MITO_EKLIPSE_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$MT_GENBANK \
					$THREADS \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

	######################################################
	##### COVERAGE STATISTICS AND PLOT FOR MT GENOME #####
	######################################################

		####################################################
		# RUN COLLECTHSMETRICS ON MT ONLY BAM FILE #########
		# USES GATK IMPLEMENTATION INSTEAD OF PICARD TOOLS #
		####################################################

			COLLECTHSMETRICS_MT ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.09-A.02-COLLECTHSMETRICS_MT"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-COLLECTHSMETRICS_MT.log" \
					-hold_jid H.09-MAKE_BAM_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.09-A.02-COLLECTHSMETRICS_MT.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$REF_GENOME \
					$MT_PICARD_INTERVAL_LIST \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

		###############################################################
		# RUN ALEX'S R SCRIPT TO GENERATE COVERAGE PLOT FOR MT GENOME #
		###############################################################

			PLOT_MT_COVERAGE ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
					$STANDARD_QUEUE_QSUB_ARG \
				-N H.09-A.02-A.01-PLOT_MT_COVERAGE"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS//$SM_TAG"-PLOT_MT_COVERAGE.log" \
					-hold_jid H.09-A.02-COLLECTHSMETRICS_MT"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/H.09-A.02-A.01_PLOT_MT_COVERAGE.sh \
					$MITO_MUTECT2_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$MT_COVERAGE_R_SCRIPT \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

###############################################################
# run steps centered on gatk's mutect2 mitochondrial workflow #
###############################################################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		CREATE_SAMPLE_ARRAY
		# run mutect2 and then filter, annotate, run haplogrep2
		MUTECT2_MT
		echo sleep 0.1s
		FILTER_MUTECT2_MT
		echo sleep 0.1s
		MASK_MUTECT2_MT
		echo sleep 0.1s
		HAPLOGREP2_MUTECT2_MT
		echo sleep 0.1s
		GNOMAD_MUTECT2_MT
		echo sleep 0.1s
		RUN_ANNOVAR_MUTECT2_MT
		echo sleep 0.1s
		FIX_ANNOVAR_MUTECT2_MT
		echo sleep 0.1s
		# run eklipse workflow
		SUBSET_BAM_MT
		echo sleep 0.1s
		RUN_EKLIPSE
		echo sleep 0.1s
		# generate coverage for mt genome
		COLLECTHSMETRICS_MT
		echo sleep 0.1s
		PLOT_MT_COVERAGE
		echo sleep 0.1s
done

#################################
##### CNV CALLING WORKFLOW ######
# USES EXOME DEPTH TO CALL CNVS #
# ANNOTATE WITH ANNOTSV #########
#################################

	########################################
	# run exomeDepth to generate CNV calls #
	########################################

		RUN_EXOME_DEPTH ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
				-N F02-RUN_EXOME_DEPTH"_"$SGE_SM_TAG"_"$PROJECT \
					-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-RUN_EXOME_DEPTH.log" \
				-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				$SCRIPT_DIR/F02_RUN_EXOME_DEPTH.sh \
					$CNV_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$SM_TAG \
					$EXOME_DEPTH_R_SCRIPT \
					$REF_PANEL_COUNTS \
					$CODING_BED \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

	################################################################
	# calculate the percent of cnv call length for each chromosome #
	################################################################

		CALCULATE_PCT_CNV_COVERAGE ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N F02-A01_PCT_CNV_COVERAGE_PER_CHR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-PCT_CNV_COVERAGE_PER_CHR.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,F02-RUN_EXOME_DEPTH"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/F02-A01_PCT_CNV_COVERAGE_PER_CHR.sh \
				$CNV_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$CODING_BED
		}

	########################################
	# run annotSV on the exomeDepth output #
	########################################

		RUN_ANNOTSV ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N F02-A02_RUN_ANNOTSV"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-RUN_ANNOTSV.log" \
			-hold_jid F02-RUN_EXOME_DEPTH"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/F02-A02_RUN_ANNOTSV.sh \
				$CNV_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	########################################################################################
	# reformat the header in the annotSV output and filter to zoom gene list if applicable #
	########################################################################################

		RUN_FORMAT_AND_ZOOM_ANNOTSV ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N F02-A02-A01_RUN_FORMAT_AND_ZOOM_ANNOTSV"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-RUN_FORMAT_AND_ZOOM_ANNOTSV.log" \
			-hold_jid F02-A02_RUN_ANNOTSV"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/F02-A02-A01_RUN_FORMAT_AND_ZOOM_ANNOTSV.sh \
				$CNV_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT \
				$ZOOM_LIST \
				$ZOOM_NAME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

##############################
# run steps for cnv workflow #
##############################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
do
	CREATE_SAMPLE_ARRAY
	RUN_EXOME_DEPTH
	echo sleep 0.1s
	CALCULATE_PCT_CNV_COVERAGE
	echo sleep 0.1s
	RUN_ANNOTSV
	echo sleep 0.1s
	RUN_FORMAT_AND_ZOOM_ANNOTSV
	echo sleep 0.1s
done

#########################################################################
##### HAPLOTYPE CALLER SCATTER ##########################################
# INPUT IS THE BAM FILE #################################################
# THE BED FILE FOR THE GVCF INTERVALS IS ################################
# THE BAIT BED PLUS THE CODING BED FILE CONCATENTATED TOGETHER ##########
# THEN A 250 BP PAD ADDED AND THEN MERGED FOR OVERLAPPING INTERVALS #####
#########################################################################

	CALL_HAPLOTYPE_CALLER ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N H.07-HAPLOTYPE_CALLER"_"$SGE_SM_TAG"_"$PROJECT"_chr"$CHROMOSOME \
			-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-HAPLOTYPE_CALLER_chr"$CHROMOSOME".log" \
		-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
		$SCRIPT_DIR/H.07_HAPLOTYPE_CALLER_SCATTER.sh \
			$GATK_3_7_0_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$SM_TAG \
			$REF_GENOME \
			$CODING_BED \
			$BAIT_BED \
			$CHROMOSOME \
			$GVCF_PAD \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

# Take the samples bait bed file and...
# create a list of unique chromosome to use as a scatter for haplotype_caller_scatter

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		|awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
do
	CREATE_SAMPLE_ARRAY
		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
			| sed -r 's/[[:space:]]+/\t/g' \
			| sed 's/chr//g' \
			| grep -v "MT" \
			| cut -f 1 \
			| sort \
			| uniq \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				collapse 1 \
			| sed 's/,/ /g');
		do
			CALL_HAPLOTYPE_CALLER
			echo sleep 0.1s
		done
done

###################################################################################################
##### HAPLOTYPE CALLER GATHER #####################################################################
# GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES AND GVCF BAM FILES INTO A SINGLE SAMPLE GVCF #
###################################################################################################

	####################################################################################################
	# create a variable to create the hold id for gathering the chromosome level gvcfs/bams per sample #
	####################################################################################################

		BUILD_HOLD_ID_PATH ()
		{
			HOLD_ID_PATH="-hold_jid "

			for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
									| sed -r 's/[[:space:]]+/\t/g' \
									| cut -f 1 \
									| sed 's/chr//g' \
									| grep -v "MT" \
									| sort \
									| uniq \
									| singularity exec $ALIGNMENT_CONTAINER datamash \
										collapse 1 \
									| sed 's/,/ /g');
				do
					HOLD_ID_PATH=$HOLD_ID_PATH"H.07-HAPLOTYPE_CALLER_"$SM_TAG"_"$PROJECT"_chr"$CHROMOSOME","
					HOLD_ID_PATH=`echo $HOLD_ID_PATH | sed 's/@/_/g'`
			done
		}

	##############################################
	# gather the per sample per chromosome gvcfs #
	##############################################

		CALL_HAPLOTYPE_CALLER_GVCF_GATHER ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.07-A.01_HAPLOTYPE_CALLER_GVCF_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG-HAPLOTYPE_CALLER_GVCF_GATHER.log \
			${HOLD_ID_PATH} \
			$SCRIPT_DIR/H.07-A.01_HAPLOTYPE_CALLER_GVCF_GATHER.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	####################################################
	# gather the per sample per chromosome hc bam file #
	####################################################

		CALL_HAPLOTYPE_CALLER_BAM_GATHER ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG-HAPLOTYPE_CALLER_BAM_GATHER.log \
			${HOLD_ID_PATH} \
			$SCRIPT_DIR/H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	########################################################
	# create a lossless HC cram, although the bam is lossy #
	########################################################

		HC_BAM_TO_CRAM ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.07-A.02-A.01_HAPLOTYPE_CALLER_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/$SM_TAG"-HC_BAM_TO_CRAM.log" \
			-hold_jid H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.07-A.02-A.01_HAPLOTYPE_CALLER_CRAM.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

##################################################################
# run steps to gather gvcfs/hc bam files and convert bam to cram #
##################################################################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		CREATE_SAMPLE_ARRAY
		BUILD_HOLD_ID_PATH
		CALL_HAPLOTYPE_CALLER_GVCF_GATHER
		echo sleep 0.1s
		CALL_HAPLOTYPE_CALLER_BAM_GATHER
		echo sleep 0.1s
		HC_BAM_TO_CRAM
		echo sleep 0.1s
done

##################################################################
##### JOINT CALLING SAMPLES IN A FAMILY WITH SET OF CONTROLS #####
##################################################################

	######################################################
	# create an array for each family/sample combination #
	######################################################

		CREATE_FAMILY_ARRAY ()
		{
			FAMILY_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} \
				$20=="'$FAMILY_ONLY'" \
				{print $1,$8,$12,$15,$16,$17,$18,$20}' \
			~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
				| sort \
				| uniq`)

				#  1  Project=the Seq Proj folder name

					PROJECT=${FAMILY_ARRAY[0]}

						################################################################################
						# 2 SKIP : FCID=flowcell that sample read group was performed on ###############
						# 3 SKIP : Lane=lane of flowcell that sample read group was performed on] ######
						# 4 SKIP : Index=sample barcode ################################################
						# 5 SKIP : Platform=type of sequencing chemistry matching SAM specification ####
						# 6 SKIP : Library_Name=library group of the sample read group #################
						# 7 SKIP : Date=should be the run set up date to match the seq run folder name #
						################################################################################

				#  8  SM_Tag=sample ID

					SM_TAG=${FAMILY_ARRAY[1]}

						SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # "@" in qsub job or holdid is not allowed

							####################################################################################
							#  9  SKIP : Center=the center/funding mechanism ###################################
							# 10  SKIP : Description=Sequencer model and/or setting (setting e.g. "Rapid-Run") #
							## Models: “HiSeq-X”,“HiSeq-4000”,“HiSeq-2500”,“HiSeq-2000”,“NextSeq-500”,“MiSeq” ##
							# 11  SKIP : Seq_Exp_ID ############################################################
							####################################################################################

				# 12  Genome_Ref=the reference genome used in the analysis pipeline

					REF_GENOME=${FAMILY_ARRAY[2]}

						########################################################
						# 13 SKIP : Operator=no standard on this, not captured #
						# 14 SKIP : Extra_VCF_Filter_Params=LEGACY, NOT USED ###
						########################################################

				# 15  TS_TV_BED_File=refseq (select) cds plus other odds and ends (.e.g. missing omim))

					TITV_BED=${FAMILY_ARRAY[3]}

				# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
				# Used for limited where to run base quality score recalibration on where to create gvcf files.

					BAIT_BED=${FAMILY_ARRAY[4]}

				# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.

					TARGET_BED=${FAMILY_ARRAY[5]}

				# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in BQSR

					DBSNP=${FAMILY_ARRAY[6]}

						#####################################################
						# 19 SKIP : KNOWN_INDEL_FILES=used for BQSR masking #
						#####################################################

				# 20 family that sample belongs to

					FAMILY=${FAMILY_ARRAY[7]}

						#######################
						# 21 SKIP : MOM #######
						# 22 SKIP : DAD #######
						# 23 SKIP : GENDER ####
						# 24 SKIP : PHENOTYPE #
						#######################

				# OLD ARRAY, DELETE LATER
					# FAMILY_PROJECT=${FAMILY_ARRAY[0]}
					# FAMILY_SGE_SAMPLE=${FAMILY_ARRAY[1]}
					# FAMILY_FAMILY=${FAMILY_ARRAY[2]}
					# FAMILY_REF_GENOME=${FAMILY_ARRAY[3]}
					# FAMILY_DBSNP=${FAMILY_ARRAY[4]}
		}

	#########################################################
	# CREATE A GVCF ".list" file for each sample per family #
	#########################################################

		CREATE_GVCF_LIST ()
		{
			awk 'BEGIN {OFS="/"} \
				$20=="'$FAMILY'" \
				{print "'$CORE_PATH'",$1,$20,$8,"GVCF",$8".g.vcf.gz"}' \
				~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
			| sort \
			| uniq \
			>| $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".gvcf.list"
		}

	############################################
	# create a list of all samples in a family #
	# *list is for gatk3, *args is for gatk4 ###
	############################################

		CREATE_FAMILY_SAMPLE_LIST ()
		{
			awk '$20=="'$FAMILY'" {print $8}' \
			~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
				| sort \
				| uniq \
			>| $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".sample.list" \
			&& cp $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".sample.list" \
			$CORE_PATH/$PROJECT/$FAMILY/$FAMILY".sample.args"
		}

	#################################################################################
	# fix common formatting problems in bed files to use for family level functions #
	# merge bait to target for gvcf creation, pad ###################################
	# create picard style interval files ############################################
	#################################################################################

		FIX_BED_FILES_FAMILY ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.10-FIX_BED_FILES"_"$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"-FIX_BED_FILES.log" \
			$SCRIPT_DIR/H.10_FIX_BED_FAMILY.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$CODING_BED \
				$TARGET_BED \
				$BAIT_BED \
				$TITV_BED \
				$CYTOBAND_BED \
				$REF_GENOME \
				$PADDING_LENGTH \
				$GVCF_PAD
		}

	###############################################################################################
	# create a hold_id variable for haplotype caller gvcf gather step for all samples in a family #
	###############################################################################################

		BUILD_HOLD_ID_PATH_GENOTYPE_GVCF ()
		{
			for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET \
				| sort \
				| uniq )
			do
				HOLD_ID_PATH="-hold_jid "
					for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
							$20=="'$FAMILY_ONLY'" \
							{print $8}' \
							~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
								| sed 's/@/_/g' \
								| sort \
								| uniq);
					do
						HOLD_ID_PATH=$HOLD_ID_PATH"H.07-A.01_HAPLOTYPE_CALLER_GVCF_GATHER_"$SAMPLE"_"$PROJECT","
					done
			done
		}

	#########################################################
	# joint calling per family per chromosome core function #
	#########################################################

		GENOTYPE_GVCF ()
			{
				echo \
				qsub \
					$QSUB_ARGS \
				-N "I.01_GENOTYPE_GVCF_SCATTER_"$FAMILY"_"$PROJECT"_chr"$CHROMOSOME \
					-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT.GENOTYPE_GVCF_chr$CHROMOSOME.log \
				$HOLD_ID_PATH \
				$SCRIPT_DIR/I.01_GENOTYPE_GVCF_SCATTER.sh \
					$GATK_3_7_0_CONTAINER \
					$CORE_PATH \
					$PROJECT \
					$FAMILY \
					$REF_GENOME \
					$DBSNP \
					$CHROMOSOME \
					$CONTROL_REPO \
					$CONTROL_DATA_SET_FILE \
					$SAMPLE_SHEET \
					$SUBMIT_STAMP
			}

	########################################
	# scatter genotype gvcfs by chromosome #
	########################################

		SCATTER_GENOTYPE_GVCF_PER_CHROMOSOME ()
		{
			for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/chr//g' \
				| grep -v "MT" \
				| cut -f 1 \
				| sort \
				| uniq \
				| singularity exec $ALIGNMENT_CONTAINER datamash \
					collapse 1 \
				| sed 's/,/ /g');
			do
				GENOTYPE_GVCF
				echo sleep 0.1s
			done
		}

#################################################################################
# run steps to do joint calling per family per set of intervals in a chromosome #
#################################################################################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
		NR>1 \
		{print $20}' \
		~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
		| sort \
		| uniq);
do
	CREATE_FAMILY_ARRAY
	CREATE_GVCF_LIST
	CREATE_FAMILY_SAMPLE_LIST
	BUILD_HOLD_ID_PATH_GENOTYPE_GVCF
	FIX_BED_FILES_FAMILY
	echo sleep 0.1s
	SCATTER_GENOTYPE_GVCF_PER_CHROMOSOME
	echo sleep 0.1s
done

########################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME GVCF FILES INTO A SINGLE FAMILY GVCF #####
########################################################################################

	########################################################################
	# create a hold_id variable for genotype gvcfs scatter step per family #
	########################################################################

		BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER ()
		{
			for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' \
				$SAMPLE_SHEET \
				| sort \
				| uniq )
			do
				HOLD_ID_PATH="-hold_jid "
				for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
					| sed -r 's/[[:space:]]+/\t/g' \
					| sed 's/chr//g' \
					| grep -v "MT" \
					| cut -f 1 \
					| sort \
					| uniq \
					| singularity exec $ALIGNMENT_CONTAINER datamash \
						collapse 1 \
					| sed 's/,/ /g');
				do
					HOLD_ID_PATH=$HOLD_ID_PATH"I.01_GENOTYPE_GVCF_SCATTER_"$FAMILY"_"$PROJECT"_chr"$CHROMOSOME","
				done
			done
		}

	###########################################################
	# gather up per chromosome genotyped vcf files per family #
	###########################################################

		CALL_GENOTYPE_GVCF_GATHER ()
		{
			echo \
			qsub \
			$QSUB_ARGS \
			-N I.01-A.01_GENOTYPE_GVCF_GATHER_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".GENOTYPE_GVCF_GATHER.log" \
			${HOLD_ID_PATH}"H.10-FIX_BED_FILES"_"$FAMILY"_"$PROJECT" \
			$SCRIPT_DIR/I.01-A.01_GENOTYPE_GVCF_GATHER.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

#####################################################
# run step to gather per chromosome per family vcfs #
#####################################################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
	~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
	| sort \
	| uniq)
do
	CREATE_FAMILY_ARRAY
	BUILD_HOLD_ID_PATH_GENOTYPE_GVCF_GATHER
	CALL_GENOTYPE_GVCF_GATHER
	echo sleep 0.1s
done

########################################################
##### DO VARIANT QUALITY SCORE RECALIBRATION ###########
# I THINK ALL OF THIS CAN BE MOVED INTO THE LOOP ABOVE #
########################################################

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		RUN_VQSR_SNP ()
		{
			echo \
			qsub \
			$QSUB_ARGS \
			-N J01_RUN_VQSR_SNP_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".RUN_VQSR_SNP.log" \
			-hold_jid I.01-A.01_GENOTYPE_GVCF_GATHER_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/J01-RUN_VARIANT_RECALIBRATOR_SNP.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$DBSNP \
				$HAPMAP \
				$OMNI_1KG \
				$HI_CONF_1KG_PHASE1_SNP \
				$SEND_TO \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	################################################
	# Run Variant Recalibrator for the INDEL model #
	################################################

		RUN_VQSR_INDEL ()
		{
			echo \
			qsub \
			$QSUB_ARGS \
			-N J02_RUN_VQSR_INDEL_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".RUN_VQSR_INDEL.log" \
			-hold_jid I.01-A.01_GENOTYPE_GVCF_GATHER_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/J02-RUN_VARIANT_RECALIBRATOR_INDEL.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$MILLS_1KG_GOLD_INDEL \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		APPLY_VQSR_SNP ()
		{
			echo \
			qsub \
			$QSUB_ARGS \
			-N K01_APPLY_VQSR_SNP_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".APPLY_VQSR_SNP.log" \
			-hold_jid J01_RUN_VQSR_SNP_$FAMILY"_"$PROJECT,J02_RUN_VQSR_INDEL_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/K01-APPLY_VARIANT_RECALIBRATION_SNP.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		APPLY_VQSR_INDEL ()
		{
			echo \
			qsub \
			$QSUB_ARGS \
			-N L01_APPLY_VQSR_INDEL_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".APPLY_VQSR_INDEL.log" \
			-hold_jid K01_APPLY_VQSR_SNP_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/L01-APPLY_VARIANT_RECALIBRATION_INDEL.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

########################
# run steps to do VQSR #
########################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
	~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
	| sort \
	| uniq)
do
	CREATE_FAMILY_ARRAY
	RUN_VQSR_SNP
	echo sleep 0.1s
	RUN_VQSR_INDEL
	echo sleep 0.1s
	APPLY_VQSR_SNP
	echo sleep 0.1s
	APPLY_VQSR_INDEL
	echo sleep 0.1s
done

################################################
##### SCATTER GATHER FOR ADDING ANNOTATION #####
################################################

	CALL_VARIANT_ANNOTATOR ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N P01_VARIANT_ANNOTATOR_$FAMILY"_"$PROJECT"_"$CHROMOSOME \
			-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".VARIANT_ANNOTATOR_$CHROMOSOME.log" \
		-hold_jid L01_APPLY_VQSR_INDEL_$FAMILY"_"$PROJECT \
		$SCRIPT_DIR/P01_VARIANT_ANNOTATOR_SCATTER.sh \
			$GATK_3_7_0_CONTAINER \
			$CORE_PATH \
			$PED_FILE \
			$PROJECT \
			$FAMILY \
			$REF_GENOME \
			$CHROMOSOME \
			$PHASE3_1KG_AUTOSOMES \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

########################
# run steps to do VQSR #
########################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
	~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
	| sort \
	| uniq);
do
	CREATE_FAMILY_ARRAY
		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
						| sed -r 's/[[:space:]]+/\t/g' \
						| sed 's/chr//g' \
						| grep -v "MT" \
						| cut -f 1 \
						| sort \
						| uniq \
						| singularity exec $ALIGNMENT_CONTAINER datamash \
							collapse 1 \
						| sed 's/,/ /g');
			do
				CALL_VARIANT_ANNOTATOR
				echo sleep 0.1s
		done
done

##############################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME ANNOTATED VCF FILES INTO A SINGLE VCF FILE #####
##############################################################################################

	######################################################
	# generate hold id from scatter of variant annotator #
	######################################################

		BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION ()
		{
			for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' \
						$SAMPLE_SHEET \
						| sort \
						| uniq )
			do
				HOLD_ID_PATH="-hold_jid "
				for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
									| sed -r 's/[[:space:]]+/\t/g' \
									| sed 's/chr//g' \
									| grep -v "MT" \
									| cut -f 1 \
									| sort \
									| uniq \
									| singularity exec $ALIGNMENT_CONTAINER datamash \
										collapse 1 \
									| sed 's/,/ /g');
				do
					HOLD_ID_PATH=$HOLD_ID_PATH"P01_VARIANT_ANNOTATOR_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
				done
			done
		}

	############################
	# variant annotator gather #
	############################

		CALL_VARIANT_ANNOTATOR_GATHER ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N P01-A01_VARIANT_ANNOTATOR_GATHER_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".MORE_VARIANT_ANNOTATOR_GATHER.log" \
			${HOLD_ID_PATH} \
			$SCRIPT_DIR/P01-A01_VARIANT_ANNOTATOR_GATHER.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#####################################################################
	# FILTER TO JUST PASSING BIALLELIC SNV SITES ON THE CODING BED FILE #
	# TEMPORARY FILE USED FOR PCA AND RELATEDNESS #######################
	#####################################################################

		CALL_PASS_BIALLELIC_SNV_COHORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N Q02-FILTER_COHORT_SNV_PASS_BIALLELIC_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".FILTER_COHORT_SNV_PASS_BIALLELIC.log" \
			-hold_jid P01-A01_VARIANT_ANNOTATOR_GATHER_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/Q02-FILTER_COHORT_SNV_PASS_BIALLELIC.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$CODING_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	################################
	# RUN PCA AND KINSHIP WORKFLOW #
	# USES KING AND PLINK ##########
	################################

		CALL_PCA_RELATEDNESS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N Q02-A01-PCA_RELATEDNESS_${FAMILY}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/${FAMILY}_${PROJECT}.PCA_RELATEDNESS.log \
			-hold_jid Q02-FILTER_COHORT_SNV_PASS_BIALLELIC_${FAMILY}_${PROJECT} \
			$SCRIPT_DIR/Q02-A01-PCA_RELATEDNESS.sh \
				$GATK_3_7_0_CONTAINER \
				$PCA_RELATEDNESS_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$PED_FILE \
				$CONTROL_PED_FILE \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

########################################################################################
########## TODO ########################################################################
########################################################################################
# ADD STEP TO FILTER ABOVE OUTPUT TO BAIT/CODING PLUS USER DEFINED PAD, PASS ONLY, ETC #
# FOR CONTROLS PLUS SAMPLES ############################################################
# THIS WOULD BE THE INPUT INTO BCFTOOLS ROH. PROBABLY ONLY A TEMP FILE #################
# MIGHT MAKE THIS WHOLE THING A SEPARATE WORKFLOW ######################################
########################################################################################

###############################################################
# run steps to do variant annotator gather and pca/relatednes #
###############################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
		~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
		| sort \
		| uniq)
	do
		CREATE_FAMILY_ARRAY
		BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION
		CALL_VARIANT_ANNOTATOR_GATHER
		echo sleep 0.1s
		CALL_PASS_BIALLELIC_SNV_COHORT
		echo sleep 0.1s
		CALL_PCA_RELATEDNESS
		echo sleep 0.1s
	done

##################################################################
##### RUNNING FILTER TO FAMILY ALL SITES BY CHROMOSOME ###########
# USE GATK4 HERE BECAUSE IT HANDLES SPANNING DELETIONS CORRECTLY #
##################################################################

	CALL_FILTER_TO_FAMILY_ALL_SITES ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N P01-A03_FILTER_TO_FAMILY_ALL_SITES_$FAMILY"_"$PROJECT"_"$CHROMOSOME \
			-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".FILTER_TO_FAMILY_ALL_SITES_$CHROMOSOME.log" \
		-hold_jid P01_VARIANT_ANNOTATOR_$FAMILY"_"$PROJECT"_"$CHROMOSOME \
		$SCRIPT_DIR/P01-A03_FILTER_TO_FAMILY_ALL_SITES_CHR.sh \
			$ALIGNMENT_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$FAMILY \
			$REF_GENOME \
			$CHROMOSOME \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

####################################################
# run steps to filter all sites vcf to family only #
####################################################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
	~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
	| sort \
	| uniq);
do
	CREATE_FAMILY_ARRAY
	for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
					| sed -r 's/[[:space:]]+/\t/g' \
					| sed 's/chr//g' \
					| grep -v "MT" \
					| cut -f 1 \
					| sort \
					| uniq \
					| singularity exec $ALIGNMENT_CONTAINER datamash \
						collapse 1 \
					| sed 's/,/ /g');
	do
		CALL_FILTER_TO_FAMILY_ALL_SITES
		echo sleep 1s
	done
done
	
#####################################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME FILTER TO FAMILY VCF FILES INTO A SINGLE VCF FILE #####
#####################################################################################################

	#################################################################################
	# create job hold id to gather up per chromosome family only all sites vcf file #
	#################################################################################

		BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF ()
		{
			for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' \
							$SAMPLE_SHEET \
							| sort \
							| uniq )
				do
					HOLD_ID_PATH="-hold_jid "
					for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
										| sed -r 's/[[:space:]]+/\t/g' \
										| sed 's/chr//g' \
										| grep -v "MT" \
										| cut -f 1 \
										| sort \
										| uniq \
										| singularity exec $ALIGNMENT_CONTAINER datamash \
											collapse 1 \
										| sed 's/,/ /g');
					do
						HOLD_ID_PATH=$HOLD_ID_PATH"P01-A03_FILTER_TO_FAMILY_ALL_SITES_"$FAMILY"_"$PROJECT"_"$CHROMOSOME","
					done
			done
		}

	######################################################
	# gather up per chromosome family only all sites vcf #
	######################################################

		CALL_FILTER_TO_FAMILY_VCF_GATHER ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N P01-A03-A01_FILTER_TO_FAMILY_ALL_SITES_GATHER_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".FILTER_TO_FAMILY_ALL_SITES_GATHER.log" \
			${HOLD_ID_PATH} \
			$SCRIPT_DIR/P01-A03-A01-FILTER_TO_FAMILY_ALL_SITES_GATHER.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	############################################################
	# filter family only all sites vcf to coding plus user pad #
	# the output for this might get moved to temp ##############
	############################################################

		CALL_FILTER_FAMILY_TO_CODING_PLUS_PAD ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N Q01_FILTER_FAMILY_CODING_PLUS_PAD_$FAMILY"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/$FAMILY"_"$PROJECT".FILTER_FAMILY_CODING_PLUS_PAD.log" \
			-hold_jid P01-A03-A01_FILTER_TO_FAMILY_ALL_SITES_GATHER_$FAMILY"_"$PROJECT \
			$SCRIPT_DIR/Q01-FILTER_FAMILY_CODING_PLUS_PAD.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$CODING_BED \
				$PADDING_LENGTH \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	############################################################
	# filter family only all sites vcf to target plus user pad #
	############################################################

		CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N Q03_FILTER_FAMILY_TARGET_PLUS_PAD_${FAMILY}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/${FAMILY}_${PROJECT}.FILTER_FAMILY_TARGET_PLUS_PAD.log \
			-hold_jid P01-A03-A01_FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY}_${PROJECT} \
			$SCRIPT_DIR/Q03-FILTER_FAMILY_TARGET_PLUS_PAD.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$TARGET_BED \
				$PADDING_LENGTH \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	############################################################
	# filter family only all sites vcf to target plus user pad #
	############################################################

		CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD_VARIANTS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N Q03-A01-FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS_${FAMILY}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/LOGS/${FAMILY}_${PROJECT}.FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS.log \
			-hold_jid Q03_FILTER_FAMILY_TARGET_PLUS_PAD_${FAMILY}_${PROJECT} \
			$SCRIPT_DIR/Q03-A01-FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

########################################################################
# run steps to gather up per chromosome family only all sites vcf file #
########################################################################

for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} {print $20}' \
	~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
	| sort \
	| uniq);
do
	CREATE_FAMILY_ARRAY
	BUILD_HOLD_ID_PATH_FILTER_TO_FAMILY_VCF
	CALL_FILTER_TO_FAMILY_VCF_GATHER
	echo sleep 0.1s
	CALL_FILTER_FAMILY_TO_CODING_PLUS_PAD
	echo sleep 0.1s
	CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD
	echo sleep 0.1s
	CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD_VARIANTS
	echo sleep 0.1s
done

#################################
### SUBSETTING TO SAMPLE VCFS ###
#################################

	#####################################################################################
	# subset sample all sites to from family coding/bait bed file plus user defined pad #
	#####################################################################################

		EXTRACT_SAMPLE_ALL_SITES ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N R01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_ALL_SITES.log \
			-hold_jid Q01_FILTER_FAMILY_CODING_PLUS_PAD_${FAMILY}_${PROJECT} \
			$SCRIPT_DIR/R01-FILTER_TO_SAMPLE_ALL_SITES.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##################################################################################
	# subset sample variant sites to from coding/bait bed file plus user defined pad #
	##################################################################################

		EXTRACT_SAMPLE_VARIANTS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N S01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_VARIANTS.log \
			-hold_jid R01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/S01-FILTER_TO_SAMPLE_VARIANTS.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#################################################################################################
	# generate vcf metrics for sample variant sites from coding/bait bed file plus user defined pad #
	#################################################################################################

		VCF_METRICS_BAIT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N T02-VCF_METRICS_BAIT_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-VCF_METRICS_BAIT.log \
			-hold_jid S01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/T02-VCF_METRICS_BAIT.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_DICT \
				$DBSNP \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#####################################################################
	# generate vcf metrics for sample variant sites from ti/tv bed file #
	#####################################################################

		VCF_METRICS_TITV ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N T03-VCF_METRICS_TITV_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-VCF_METRICS_TITV.log \
			-hold_jid S01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/T03-VCF_METRICS_TITV.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_DICT \
				$TITV_BED \
				$DBSNP_129 \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#########################################################################
	# subset sample to all sites from target bed file plus user defined pad #
	#########################################################################

		EXTRACT_SAMPLE_ALL_SITES_ON_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N T01-FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_ALL_SITES_TARGET.log \
			-hold_jid R01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/T01-FILTER_TO_SAMPLE_ALL_SITES_TARGET.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$TARGET_BED \
				$PADDING_LENGTH \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	########################################################################
	# subset sample variant sites to target bed file plus user defined pad #
	########################################################################

		EXTRACT_SAMPLE_VARIANTS_ON_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N U01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_VARIANTS_TARGET.log \
			-hold_jid T01-FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/U01-FILTER_TO_SAMPLE_VARIANTS_TARGET.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#####################################################################################
	# decompose on target plus user defined pad sample variant sites to target bed file #
	#####################################################################################

		DECOMPOSE_SAMPLE_VARIANTS_ON_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N U01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-DECOMPOSE_SAMPLE_VARIANTS_TARGET.log \
			-hold_jid U01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/U01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET.sh \
				$VT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#############################################################################
	# run annovar on decomposed on target plus user defined pad sample vcf file #
	#############################################################################

		RUN_ANNOVAR_ON_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N U01-A01-A01-RUN_ANNOVAR_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-RUN_ANNOVAR_TARGET.log \
			-hold_jid U01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/U01-A01-A01-RUN_ANNOVAR_TARGET.sh \
				$ANNOVAR_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$ANNOVAR_DATABASE_FILE \
				$ANNOVAR_REF_BUILD \
				$ANNOVAR_INFO_FIELD_KEYS \
				$ANNOVAR_HEADER_MAPPINGS \
				$ANNOVAR_VCF_COLUMNS \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	############################################################################################
	# generate vcf metrics for sample variant sites from target bed file plus user defined pad #
	############################################################################################

		VCF_METRICS_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N V01-VCF_METRICS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-VCF_METRICS_TARGET.log \
			-hold_jid U01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			$SCRIPT_DIR/V01-VCF_METRICS_TARGET.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FAMILY \
				$SM_TAG \
				$REF_DICT \
				$DBSNP \
				$THREADS \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##################################
	# QC REPORT PREP FOR EACH SAMPLE #
	##################################

QC_REPORT_PREP ()
{
echo \
qsub \
	$QSUB_ARGS \
-N X01-QC_REPORT_PREP_${SGE_SM_TAG}_${PROJECT} \
	-o $CORE_PATH/$PROJECT/$FAMILY/$SM_TAG/LOGS/${SM_TAG}-QC_REPORT_PREP.log \
-hold_jid \
V01-VCF_METRICS_TARGET_${SGE_SM_TAG}_${PROJECT},\
T03-VCF_METRICS_TITV_${SGE_SM_TAG}_${PROJECT},\
T02-VCF_METRICS_BAIT_${SGE_SM_TAG}_${PROJECT},\
H.05-A.01_CHROM_DEPTH_${SGE_SM_TAG}_${PROJECT},\
H.04-A.01-RUN_VERIFYBAMID_${SGE_SM_TAG}_${PROJECT},\
H.02-COLLECT_HS_METRICS_${SGE_SM_TAG}_${PROJECT},\
H.01-COLLECT_MULTIPLE_METRICS_${SGE_SM_TAG}_${PROJECT},\
F.01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
$SCRIPT_DIR/X01-QC_REPORT_PREP.sh \
	$ALIGNMENT_CONTAINER \
	$CORE_PATH \
	$PROJECT \
	$FAMILY \
	$SM_TAG \
	$FATHER \
	$MOTHER \
	$GENDER \
	$PHENOTYPE
}

###############################
# run vcf sample subset steps #
###############################

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
	| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
	| awk 'BEGIN {FS=","} NR>1 {print $8}' \
	| sort \
	| uniq );
do
	CREATE_SAMPLE_ARRAY
	EXTRACT_SAMPLE_ALL_SITES
	echo sleep 0.1s
	EXTRACT_SAMPLE_VARIANTS
	echo sleep 0.1s
	VCF_METRICS_BAIT
	echo sleep 0.1s
	VCF_METRICS_TITV
	echo sleep 0.1s
	EXTRACT_SAMPLE_ALL_SITES_ON_TARGET
	echo sleep 0.1s
	EXTRACT_SAMPLE_VARIANTS_ON_TARGET
	echo sleep 0.1s
	DECOMPOSE_SAMPLE_VARIANTS_ON_TARGET
	echo sleep 0.1s
	RUN_ANNOVAR_ON_TARGET
	echo sleep 0.1s
	VCF_METRICS_TARGET
	echo sleep 0.1s
	QC_REPORT_PREP
	echo sleep 0.1s
done

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

