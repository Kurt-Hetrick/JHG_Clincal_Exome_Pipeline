#!/usr/bin/env bash

###################
# INPUT VARIABLES #
###################

	SAMPLE_SHEET=$1
	PED_FILE=$2
	PADDING_LENGTH=$3 # optional. if no 3rd argument present then the default is 10
	# THIS PAD IS FOR SLICING

		if [[ ! ${PADDING_LENGTH} ]]
			then
			PADDING_LENGTH="10"
		fi

	QUEUE_LIST=$4 # optional. if no 4th argument present then the default is cgc.q
		# if you want to set this then you need to set the 3rd argument as well (even to the default)

		if [[ ! ${QUEUE_LIST} ]]
			then
			QUEUE_LIST="cgc.q"
		fi

	PRIORITY=$5 # optional. if no 5th argument present then the default is -15.
		# if you want to set this then you need to set the 3rd and 4th argument as well (even to the default)

			if [[ ! ${PRIORITY} ]]
				then
				PRIORITY="-15"
			fi

	THREADS=$6 # optional. if no 6th argument present then default is 6.
		# if you want to set this then you need to set 3rd,4th and 5th argument as well (even to default)

			if [[ ! ${THREADS} ]]
				then
				THREADS="6"
			fi

########################################################################
# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED #
########################################################################

	SUBMITTER_SCRIPT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

	SCRIPT_DIR=${SUBMITTER_SCRIPT_PATH}/scripts

##################
# CORE VARIABLES #
##################

	# GVCF PAD. CURRENTLY KEEPING THIS AS A STATIC VARIABLE

		GVCF_PAD="250"

	## This will always put the current working directory in front of any directory for PATH
	## added /bin for RHEL6

		export PATH=".:${PATH}:/bin"

	# where the input/output sequencing data will be located.

		CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"

	# Directory where NovaSeqa runs are located.

		NOVASEQ_REPO="/mnt/instrument_files/novaseq"

	# used for tracking in the read group header of the cram file

		PIPELINE_VERSION=$(git --git-dir=${SCRIPT_DIR}/../.git --work-tree=${SCRIPT_DIR}/.. log --pretty=format:'%h' -n 1)

	# load gcc for programs like verifyBamID
	## this will get pushed out to all of the compute nodes since I specify env var to pushed out with qsub

		module load gcc/7.2.0

	# explicitly setting this b/c not everybody has had the $HOME directory transferred
	# and I'm not going to through and figure out who does and does not have this set correctly

		umask 0007

	# SUBMIT TIMESTAMP

		SUBMIT_STAMP=$(date '+%s')

	# SUBMITTER_ID

		SUBMITTER_ID=$(whoami)

	# grab submitter's name

		PERSON_NAME=$(getent passwd | awk 'BEGIN {FS=":"} $1=="'${SUBMITTER_ID}'" {print $5}')

	# grab email addy

		SEND_TO=$(cat ${SCRIPT_DIR}/../email_lists.txt)

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
				QSUB_ARGS=${QSUB_ARGS}" -cwd" \
				QSUB_ARGS=${QSUB_ARGS}" -V" \
				QSUB_ARGS=${QSUB_ARGS}" -v SINGULARITY_BINDPATH=/mnt:/mnt" \
				QSUB_ARGS=${QSUB_ARGS}" -q ${QUEUE_LIST}" \
				QSUB_ARGS=${QSUB_ARGS}" -p ${PRIORITY}" \
				QSUB_ARGS=${QSUB_ARGS}" -j y"

			# qsub args for magick package in R (imgmagick_merge.r)
			# image packages will use all cpu threads by default.
			# to configure set env variable to desired thread count.

				IMGMAGICK_QSUB_ARGS=${QSUB_ARGS}" -v MAGICK_THREAD_LIMIT=${THREADS}"

#####################
# PIPELINE PROGRAMS #
#####################

	############################
	# BASE PIPELINE CONTAINERS #
	############################

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

	################################################################
	# MITOCHONDRIA ANALYSIS CONTAINERS AND AUXILIARY SCRIPTS/FILES #
	################################################################

		MITO_MUTECT2_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/mito_mutect2-4.1.3.0.0.simg"
			# uses broadinstitute/gatk:4.1.3.0 as the base image (as /gatk/gatk.jar)
				# added
					# bcftools-1.10.2
					# haplogrep-2.1.20.jar (as /jars/haplogrep-2.1.20.jar)
					# annovar

		MITO_EKLIPSE_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/mito_eklipse-master-c25931b.0.simg"
			# https://github.com/dooguypapua/eKLIPse AND all of its dependencies

		MITO_MAGICK_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/mito_magick-6.8.9.9.0.simg"
			# magick package for R. see dockerfile for details.

		EKLIPSE_CIRCOS_LEGEND="${SCRIPT_DIR}/circos_legend.png"

		EKLIPSE_FORMAT_CIRCOS_PLOT_R_SCRIPT="${SCRIPT_DIR}/imgmagick_merge.r"

		MT_COVERAGE_R_SCRIPT="${SCRIPT_DIR}/mito_coverage_graph.r"

	#################################################
	# CNV ANALYSIS CONTAINERS AND AUXILIARY SCRIPTS #
	#################################################

		CNV_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/cnv_exomedepth-dev.simg"

		EXOME_DEPTH_R_SCRIPT="${SCRIPT_DIR}/runExomeDepth.r"

		FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT="${SCRIPT_DIR}/FORMAT_AND_ZOOM_ANNOTSV.r"

	#################################
	# PCA AND RELATEDNESS CONTAINER #
	#################################

		PCA_RELATEDNESS_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/pca-relatedness-0.0.1.simg"

	#################################################
	# WES ANNOVAR CONTAINERS, PARAMETERS AND INPUTS #
	#################################################

		VT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/vt-0.5772.ca352e2c.0.simg"

		ANNOVAR_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/containers/annovarwrangler-20210126.simg"

		# ANNOVAR PARAMETERS AND INPUTS

			ANNOVAR_DATABASE_FILE="${SCRIPT_DIR}/../resources/CFTR.final.csv"
			ANNOVAR_REF_BUILD="hg19"

			ANNOVAR_INFO_FIELD_KEYS="VariantType," \
				ANNOVAR_INFO_FIELD_KEYS=${ANNOVAR_INFO_FIELD_KEYS}"DP" \

			ANNOVAR_HEADER_MAPPINGS="af=gnomad211_exome_AF," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_popmax=gnomad211_exome_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_male=gnomad211_exome_AF_male," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_female=gnomad211_exome_AF_female," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_raw=gnomad211_exome_AF_raw," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_afr=gnomad211_exome_AF_afr," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_sas=gnomad211_exome_AF_sas," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_amr=gnomad211_exome_AF_amr," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_eas=gnomad211_exome_AF_eas," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_nfe=gnomad211_exome_AF_nfe," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_fin=gnomad211_exome_AF_fin," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_asj=gnomad211_exome_AF_asj," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"af_oth=gnomad211_exome_AF_oth," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_topmed_af_popmax=gnomad211_exome_non_topmed_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_neuro_af_popmax=gnomad211_exome_non_neuro_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_cancer_af_popmax=gnomad211_exome_non_cancer_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"controls_af_popmax=gnomad211_exome_controls_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF=gnomad211_genome_AF," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_popmax=gnomad211_genome_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_male=gnomad211_genome_AF_male," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_female=gnomad211_genome_AF_female," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_raw=gnomad211_genome_AF_raw," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_afr=gnomad211_genome_AF_afr," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_sas=gnomad211_genome_AF_sas," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_amr=gnomad211_genome_AF_amr," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_eas=gnomad211_genome_AF_eas," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_nfe=gnomad211_genome_AF_nfe," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_fin=gnomad211_genome_AF_fin," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_asj=gnomad211_genome_AF_asj," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"AF_oth=gnomad211_genome_AF_oth," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_topmed_AF_popmax=gnomad211_genome_non_topmed_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_neuro_AF_popmax=gnomad211_genome_non_neuro_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"non_cancer_AF_popmax=gnomad211_genome_non_cancer_AF_popmax," \
				ANNOVAR_HEADER_MAPPINGS=${ANNOVAR_HEADER_MAPPINGS}"controls_AF_popmax=gnomad211_genome_controls_AF_popmax"

				ANNOVAR_VCF_COLUMNS="CHROM,"
					ANNOVAR_VCF_COLUMNS=${ANNOVAR_VCF_COLUMNS}"POS,"
					ANNOVAR_VCF_COLUMNS=${ANNOVAR_VCF_COLUMNS}"REF,"
					ANNOVAR_VCF_COLUMNS=${ANNOVAR_VCF_COLUMNS}"ALT"

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
		UCSC_REPEATMASK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINES/JHG_Clinical_Exome_Pipeline_Phase2/resources/ucsc_grch37_repeatmasker.sorted_no_alt_MT.bed"
			# sortBed -i ucsc_grch37_repeatmasker.bed \
			# | awk '$1!~"_"&&$1!~"chrM"' \
			# | sed 's/^chr//g' \
			# > ucsc_grch37_repeatmasker.sorted_no_alt_MT.bed
		MDUST_REPEATMASK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINES/JHG_Clinical_Exome_Pipeline_Phase2/resources/LCR-hs37d5.bed"
			# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271055/
			# https://github.com/lh3/varcmp/tree/master/scripts

	# where the control data set resides.

		CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Data/TWIST_CONTROL_SET1.200601_PIPELINE_2_0_0"
		CONTROL_PED_FILE="${CONTROL_REPO}/TWIST_CONTROL_SET1.200601.ped"

	# SFAFASFA

		CONTROL_DATA_SET_FILE="CGC_CONTROL_SET_3_7.g.vcf.gz"

	# Mitochondrial Pipline

		MT_PICARD_INTERVAL_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/MT.interval_list"
		MT_MASK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/hg37_MT_blacklist_sites.hg37.MT.bed"
		GNOMAD_MT="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/GRCh37_MT_gnomAD.vcf.gz"
		# ANNOVAR_MT_DB_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/annovar_db/"
		ANNOVAR_MT_DB_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/annovar_db/2021_02_02/annovar/humandb"
		MT_GENBANK="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/MITO/NC_012920.1.gb"

	# CNV calling workflow

		## REF_PANEL_COUNTS USED IN EXOME DEPTH IS SEX SPECIFIC.
		## DETERMINED WHEN PARSING GENDER FROM PED FILE DURING CREATE_SAMPLE_ARRAY
		## THE THREE RDA FILES BELOW GET REASSIGNED TO ${REF_PANEL_COUNTS} depending on what the gender is.

			# read count from female reference panel, won't need to change unless changes in bed file or reference samples

				REF_PANEL_FEMALE_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountFemaleUniqBed22.rda"

			# read count from male reference panel, won't need to change unless changes in bed file or reference samples

				REF_PANEL_MALE_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountMaleUniqBed26.rda"

			# if subject sex is not specified as 'm' or 'f', it will use count of all sample

				REF_PANEL_ALL_READ_COUNT_RDA="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/CNV/refCountAllUniqBed48.rda"

############################################################################
##### PIPELINE AND PROJECT SET-UP ##########################################
############################################################################
##### MERGE SAMPLE_SHEET AND PED FILE AND CREATE A SAMPLE LEVEL ARRAY ######
##### ARRAY IS USED TO PASS VARIABLES FOR SAMPLE LEVEL PROCESSES ###########
##### MAKE A DIRECTORY TREE ################################################
##### FIX BED FILES USED FOR EACH SAMPLE ###################################
##### FIX BED FILES USED FOR EACH FAMILY ###################################
##### CREATE LISTS FOR SAMPLES IN A FAMILY (USED FOR JOINT CALLING, ETC) ###
############################################################################

########################################
### MERGE SAMPLE SHEET WITH PED FILE ###
########################################

	# make a directory in user home directory

		mkdir -p ~/CGC_PIPELINE_TEMP

	# create variables using the base name for the sample sheet and ped file

		MANIFEST_PREFIX=$(basename ${SAMPLE_SHEET} .csv)
		PED_PREFIX=$(basename ${PED_FILE} .ped)

	# fix any commonly seen formatting issues in the sample sheet

		FORMAT_MANIFEST ()
		{
			awk 1 ${SAMPLE_SHEET} \
				| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
				| awk 'NR>1' \
				| sed 's/,/\t/g' \
				| sort -k 8,8 \
			>| ~/CGC_PIPELINE_TEMP/SORTED.${MANIFEST_PREFIX}.txt
		}

	# merge the sample sheet with the ped file

		MERGE_PED_MANIFEST ()
		{
			awk 1 ${PED_FILE} \
				| sed 's/\r//g' \
				| sort -k 2,2 \
				| join -1 8 -2 2 -e '-'  -t $'\t' \
				-o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
			~/CGC_PIPELINE_TEMP/SORTED.${MANIFEST_PREFIX}.txt /dev/stdin \
			>| ~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt
		}

#############################################################################
### CREATE_SAMPLE_ARRAY #####################################################
### create an array from values of the merged sample sheet and ped file #####
### set ${REF_PANEL_COUNTS} based on gender #################################
#############################################################################

	CREATE_SAMPLE_ARRAY ()
	{
		SAMPLE_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} \
			$8=="'${SAMPLE}'" \
			{split($19,INDEL,";"); \
			print $1,$8,$9,$10,$11,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],\
			$20,$21,$22,$23,$24}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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

					# If there is an @ in the qsub or holdId name it breaks

						SGE_SM_TAG=$(echo ${SM_TAG} | sed 's/@/_/g')

			#  9  Center=the center/funding mechanism

				CENTER=${SAMPLE_ARRAY[2]}

			# 10  Description=Sequencer model and/or setting (setting e.g. "Rapid-Run")
			## Models: “HiSeq-X”,“HiSeq-4000”,“HiSeq-2500”,“HiSeq-2000”,“NextSeq-500”,“MiSeq”

				SEQUENCER_MODEL=${SAMPLE_ARRAY[3]}

			# 11  Seq_Exp_ID=Zoom Gene List for filtering ExomeDepth output by gene symbol

				ZOOM_LIST=${SAMPLE_ARRAY[4]}

					# if the zoom list file exists than the output file prefix is the input file prefix before .GeneList

						if [ -f ${ZOOM_LIST} ]
							then ZOOM_NAME=$(basename ${ZOOM_LIST} | sed 's/.GeneList.[0-9]*.csv//g')
							else ZOOM_NAME="NA"
						fi

			# 12  Genome_Ref=the reference genome used in the analysis pipeline

				REF_GENOME=${SAMPLE_ARRAY[5]}

					# REFERENCE DICTIONARY IS A SUMMARY OF EACH CONTIG. PAIRED WITH REF GENOME

						REF_DICT=$(echo ${REF_GENOME} | sed 's/fasta$/dict/g; s/fa$/dict/g')

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

				# set ${REF_PANEL_COUNTS} USED IN EXOMEDEPTH TO THE SEX SPECIFIC ONE

					if [[ ${GENDER} = "1" ]];
						then REF_PANEL_COUNTS=${REF_PANEL_MALE_READ_COUNT_RDA}
						elif [[ ${GENDER} = "2" ]];
							then REF_PANEL_COUNTS=${REF_PANEL_FEMALE_READ_COUNT_RDA}
						else
							REF_PANEL_COUNTS=${REF_PANEL_ALL_READ_COUNT_RDA}
					fi

			# 24 PHENOTYPE

				PHENOTYPE=${SAMPLE_ARRAY[16]}
	}

######################################
### PROJECT DIRECTORY TREE CREATOR ###
######################################

	MAKE_PROJ_DIR_TREE ()
	{
		mkdir -p \
		${CORE_PATH}/${PROJECT}/{FASTQ,LOGS,COMMAND_LINES,REPORTS} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/{LOGS,PCA,RELATEDNESS,ROH} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/VCF/{RAW,VQSR} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/{CNV_OUTPUT,CRAM,GVCF,HC_CRAM,LOGS} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/MT_OUTPUT/{COLLECTHSMETRICS_MT,MUTECT2_MT,HAPLOGROUPS,ANNOVAR_MT,EKLIPSE} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/{ALIGNMENT_SUMMARY,ANEUPLOIDY_CHECK,ANNOVAR,ERROR_SUMMARY,PICARD_DUPLICATES,QC_REPORT_PREP,QUALITY_YIELD,RG_HEADER,TI_TV,VCF_METRICS,VERIFYBAMID,VERIFYBAMID_AUTO} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/{TARGET_PADDED,CODING_PADDED} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/INSERT_SIZE/{METRICS,PDF} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
		${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		${CORE_PATH}/${PROJECT}/TEMP/{KING,PLINK,VCF_PREP} \
		${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_ANNOVAR_MT \
		${CORE_PATH}/${PROJECT}/TEMP/${SM_TAG}_ANNOVAR_TARGET
	}

############################################################################
### combine above functions into one...this is probably not necessary... ###
############################################################################

	SETUP_PROJECT ()
	{
		FORMAT_MANIFEST
		MERGE_PED_MANIFEST
		CREATE_SAMPLE_ARRAY
		MAKE_PROJ_DIR_TREE
		echo Project started at `date` >| ${CORE_PATH}/${PROJECT}/REPORTS/PROJECT_START_END_TIMESTAMP.txt
	}

###################################################
### fix common formatting problems in bed files ###
### merge bait to target for gvcf creation, pad ###
### create picard style interval files ############
### DO PER SAMPLE #################################
###################################################

	FIX_BED_FILES_SAMPLE ()
	{
		echo \
		qsub \
			${QSUB_ARGS} \
		-N A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FIX_BED_FILES.log \
		${SCRIPT_DIR}/A02-FIX_BED_FILES_SAMPLE.sh \
			${ALIGNMENT_CONTAINER} \
			${CORE_PATH} \
			${PROJECT} \
			${SM_TAG} \
			${CODING_BED} \
			${TARGET_BED} \
			${BAIT_BED} \
			${TITV_BED} \
			${CYTOBAND_BED} \
			${REF_GENOME} \
			${REF_DICT} \
			${PADDING_LENGTH} \
			${GVCF_PAD}
	}

######################################################
### CREATE_FAMILY_ARRAY ##############################
# create an array for each family/sample combination #
######################################################

	CREATE_FAMILY_ARRAY ()
	{
		FAMILY_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} \
			$20=="'${FAMILY_ONLY}'" \
			{print $1,$8,$12,$15,$16,$17,$18,$20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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

					# "@" in qsub job or holdid is not allowed

						SGE_SM_TAG=$(echo ${SM_TAG} | sed 's/@/_/g')

							####################################################################################
							#  9  SKIP : Center=the center/funding mechanism ###################################
							# 10  SKIP : Description=Sequencer model and/or setting (setting e.g. "Rapid-Run") #
							## Models: “HiSeq-X”,“HiSeq-4000”,“HiSeq-2500”,“HiSeq-2000”,“NextSeq-500”,“MiSeq” ##
							# 11  SKIP : Seq_Exp_ID ############################################################
							####################################################################################

			# 12  Genome_Ref=the reference genome used in the analysis pipeline

				REF_GENOME=${FAMILY_ARRAY[2]}

					# REFERENCE DICTIONARY IS A SUMMARY OF EACH CONTIG. PAIRED WITH REF GENOME

						REF_DICT=$(echo ${REF_GENOME} | sed 's/fasta$/dict/g; s/fa$/dict/g')

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
	}

#############################################################
### CREATE A GVCF ".list" file for each sample per family ###
#############################################################

	CREATE_GVCF_LIST ()
	{
		awk 'BEGIN {FS="\t"; OFS="/"} \
			$20=="'${FAMILY}'" \
			{print "'${CORE_PATH}'",$1,$20,$8,"GVCF",$8".g.vcf.gz"}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq \
		>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${FAMILY}.gvcf.list
	}

################################################
### create a list of all samples in a family ###
### *list is for gatk3, *args is for gatk4 #####
################################################

	CREATE_FAMILY_SAMPLE_LIST ()
	{
		awk 'BEGIN {FS="\t"; OFS="\t"} \
			$20=="'${FAMILY}'" \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq \
		>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${FAMILY}.sample.list \
			&& cp ${CORE_PATH}/${PROJECT}/${FAMILY}/${FAMILY}.sample.list \
			${CORE_PATH}/${PROJECT}/${FAMILY}/${FAMILY}.sample.args
	}

#####################################################################################
### fix common formatting problems in bed files to use for family level functions ###
### merge bait to target for gvcf creation, pad #####################################
### create picard style interval files ##############################################
#####################################################################################

	FIX_BED_FILES_FAMILY ()
	{
		echo \
		qsub \
			${QSUB_ARGS} \
		-N A03-FIX_BED_FILES_${FAMILY}_${PROJECT} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}-FIX_BED_FILES.log \
		${SCRIPT_DIR}/A03-FIX_BED_FILES_FAMILY.sh \
			${ALIGNMENT_CONTAINER} \
			${CORE_PATH} \
			${PROJECT} \
			${FAMILY} \
			${CODING_BED} \
			${TARGET_BED} \
			${BAIT_BED} \
			${TITV_BED} \
			${CYTOBAND_BED} \
			${REF_GENOME} \
			${REF_DICT} \
			${PADDING_LENGTH} \
			${GVCF_PAD}
	}

############################################
# RUN STEPS FOR PIPELINE AND PROJECT SETUP #
############################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		SETUP_PROJECT
		FIX_BED_FILES_SAMPLE
		echo sleep 0.1s
	done

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_FAMILY_ARRAY
		CREATE_GVCF_LIST
		CREATE_FAMILY_SAMPLE_LIST
		FIX_BED_FILES_FAMILY
		echo sleep 0.1s
	done

#######################################################################################
##### CRAM FILE GENERATION ############################################################
#######################################################################################
# NOTE: THE CRAM FILE IS THE END PRODUCT BUT THE BAM FILE IS USED FOR OTHER PROCESSES #
# SOME PROGRAMS CAN'T TAKE IN CRAM AS AN INPUT ########################################
#######################################################################################

#################################################################################################
### CREATE_PLATFROM_UNIT_ARRAY ##################################################################
### create an array at the platform unit level so that bwa mem can add metadata to the header ###
#################################################################################################

	CREATE_PLATFORM_UNIT_ARRAY ()
	{
		PLATFORM_UNIT_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} \
			$8$2$3$4=="'${PLATFORM_UNIT}'" \
			{split($19,INDEL,";"); \
			print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],\
			$20,$21,$22,$23,$24}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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

						SGE_SM_TAG=$(echo ${SM_TAG} | sed 's/@/_/g')

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
### RUN_BWA ############################################################
### Use bwa mem to do the alignments; ##################################
### pipe to samblaster to add mate tags; ###############################
### pipe to picard's AddOrReplaceReadGroups to handle the bam header ###
########################################################################

	RUN_BWA ()
	{
		echo \
		qsub \
			${QSUB_ARGS} \
		-N A01-BWA_${SGE_SM_TAG}_${FCID}_${LANE}_${INDEX} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}_${FCID}_${LANE}_${INDEX}-BWA.log \
		${SCRIPT_DIR}/A01-BWA.sh \
			${ALIGNMENT_CONTAINER} \
			${CORE_PATH} \
			${PROJECT} \
			${FCID} \
			${LANE} \
			${INDEX} \
			${PLATFORM} \
			${LIBRARY} \
			${RUN_DATE} \
			${SM_TAG} \
			${CENTER} \
			${SEQUENCER_MODEL} \
			${REF_GENOME} \
			${PIPELINE_VERSION} \
			${BAIT_BED} \
			${TARGET_BED} \
			${TITV_BED} \
			${NOVASEQ_REPO} \
			${THREADS} \
			${SAMPLE_SHEET} \
			${SUBMIT_STAMP}
	}

#############################
# RUN STEPS TO RUN BWA, ETC #
#############################

	for PLATFORM_UNIT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8$2$3$4}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq );
	do
		CREATE_PLATFORM_UNIT_ARRAY
		RUN_BWA
		echo sleep 0.1s
	done

#########################################################################################
### MARK_DUPLICATES #####################################################################
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

			awk 'BEGIN {FS="\t"; OFS="\t"} \
				{print $1,$20,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam",$8,$10}' \
			~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| awk 'BEGIN {OFS="\t"} \
				{sub(/@/,"_",$6)} \
				{print $1,$2,$3,$4,$5,$6,$7}' \
			| sort -k 1,1 -k 2,2 -k 3,3 -k 4,4 -k 7,7 \
			| uniq \
			| singularity exec ${ALIGNMENT_CONTAINER} datamash \
				-s \
				-g 1,2,3 \
				collapse 4 \
				collapse 5 \
				unique 6 \
				unique 7 \
			| awk 'BEGIN {FS="\t"} \
				gsub(/,/,",A01-BWA_"$6"_",$4) \
				gsub(/,/,",INPUT=" "'${CORE_PATH}'" "/" $1"/TEMP/",$5) \
				{print "qsub",\
				"-S /bin/bash",\
				"-cwd",\
				"-V",\
				"-v SINGULARITY_BINDPATH=/mnt:/mnt",\
				"-q","'${QUEUE_LIST}'",\
				"-p","'${PRIORITY}'",\
				"-j y",\
				"-N","B01-MARK_DUPLICATES_"$6"_"$1,\
				"-o","'${CORE_PATH}'/"$1"/"$2"/"$3"/LOGS/"$3"-MARK_DUPLICATES.log",\
				"-hold_jid","A01-BWA_"$6"_"$4, \
				"'${SCRIPT_DIR}'""/B01-MARK_DUPLICATES.sh",\
				"'${ALIGNMENT_CONTAINER}'",\
				"'${CORE_PATH}'",\
				$1,\
				$2,\
				$3,\
				$7,\
				"'${THREADS}'",\
				"'${SAMPLE_SHEET}'",\
				"'${SUBMIT_STAMP}'",\
				"INPUT=" "'${CORE_PATH}'" "/" $1"/TEMP/"$5"\n""sleep 0.1s"}'

########################################
### FINISH UP CREATING THE CRAM FILE ###
# bqsr then convert cram to BAM ########
########################################

	#######################################
	# run bqsr on the using bait bed file #
	#######################################

		PERFORM_BQSR ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N C01-PERFORM_BQSR_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-PERFORM_BQSR.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},B01-MARK_DUPLICATES_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/C01-PERFORM_BQSR.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${KNOWN_INDEL_1} \
				${KNOWN_INDEL_2} \
				${DBSNP} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
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
				${QSUB_ARGS} \
			-N D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-APPLY_BQSR.log \
			-hold_jid C01-PERFORM_BQSR_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/D01-APPLY_BQSR.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#####################################################
	# create a lossless cram, although the bam is lossy #
	#####################################################

		BAM_TO_CRAM ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-BAM_TO_CRAM.log \
			-hold_jid D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E01-BAM_TO_CRAM.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

######################################
# RUN STEPS FOR CRAM FILE GENERATION #
######################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_SAMPLE_ARRAY
		PERFORM_BQSR
		echo sleep 0.1s
		APPLY_BQSR
		echo sleep 0.1s
		BAM_TO_CRAM
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
			${QSUB_ARGS} \
		-N E02-HAPLOTYPE_CALLER_${SGE_SM_TAG}_${PROJECT}_chr${CHROMOSOME} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-HAPLOTYPE_CALLER_chr${CHROMOSOME}.log \
		-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
		${SCRIPT_DIR}/E02-HAPLOTYPE_CALLER_SCATTER.sh \
			${GATK_3_7_0_CONTAINER} \
			${CORE_PATH} \
			${PROJECT} \
			${SM_TAG} \
			${REF_GENOME} \
			${CODING_BED} \
			${BAIT_BED} \
			${CHROMOSOME} \
			${GVCF_PAD} \
			${SAMPLE_SHEET} \
			${SUBMIT_STAMP}
	}

#######################################################################################
# RUN STEPS FOR HAPLOTYPE CALLER SCATTER ##############################################
# Take the samples bait bed file and ##################################################
# create a list of unique chromosome to use as a scatter for haplotype_caller_scatter #
#######################################################################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_SAMPLE_ARRAY
			for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/chr//g' \
				| egrep "^[0-9]|^X|^Y" \
				| cut -f 1 \
				| sort -V \
				| uniq \
				| singularity exec ${ALIGNMENT_CONTAINER} datamash \
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

			for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/chr//g' \
				| egrep "^[0-9]|^X|^Y" \
				| cut -f 1 \
				| sort -V \
				| uniq \
				| singularity exec ${ALIGNMENT_CONTAINER} datamash \
					collapse 1 \
				| sed 's/,/ /g');
			do
				HOLD_ID_PATH="${HOLD_ID_PATH}E02-HAPLOTYPE_CALLER_${SM_TAG}_${PROJECT}_chr${CHROMOSOME},"

				HOLD_ID_PATH=`echo ${HOLD_ID_PATH} | sed 's/@/_/g'`
			done
		}

	##############################################
	# gather the per sample per chromosome gvcfs #
	##############################################

		CALL_HAPLOTYPE_CALLER_GVCF_GATHER ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E02-A01-HAPLOTYPE_CALLER_GVCF_GATHER_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-HAPLOTYPE_CALLER_GVCF_GATHER.log \
			${HOLD_ID_PATH} \
			${SCRIPT_DIR}/E02-A01-HAPLOTYPE_CALLER_GVCF_GATHER.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	####################################################
	# gather the per sample per chromosome hc bam file #
	####################################################

		CALL_HAPLOTYPE_CALLER_BAM_GATHER ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E02-A02-HAPLOTYPE_CALLER_BAM_GATHER_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-HAPLOTYPE_CALLER_BAM_GATHER.log \
			${HOLD_ID_PATH} \
			${SCRIPT_DIR}/E02-A02-HAPLOTYPE_CALLER_BAM_GATHER.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${SM_TAG} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	########################################################
	# create a lossless HC cram, although the bam is lossy #
	########################################################

		HC_BAM_TO_CRAM ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E02-A02-A01-HAPLOTYPE_CALLER_CRAM_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-HC_BAM_TO_CRAM.log \
			-hold_jid E02-A02-HAPLOTYPE_CALLER_BAM_GATHER_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E02-A02-A01-HAPLOTYPE_CALLER_CRAM.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

##################################################################
# RUN STEPS TO GATHER GVCFS/HC BAM FILES AND CONVERT BAM TO CRAM #
##################################################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
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

#############################################
##### MITOCHONDRIAL WORKFLOW ################
#############################################
# RUN MUTECT2 TO CALL SNVS AND SMALL INDELS #
# GENERATE COVERAGE STATS ###################
# RUN EKLIPSE FOR LARGER DELETIONS ##########
#############################################

	#####################################
	### MUTECT2 IN MITO MODE WORKFLOW ###
	### WORKS ON FULL BAM FILE ##########
	#####################################

		#####################################################
		# run mutect2 in mitochondria mode on full bam file #
		# this runs MUCH slower on non-avx machines #########
		#####################################################

			MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-MUTECT2_MT.log \
				-hold_jid D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-MUTECT2_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${REF_GENOME} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		#######################################
		# apply filters to mutect2 vcf output #
		#######################################

			FILTER_MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A01-FILTER_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_MUTECT2_MT.log \
				-hold_jid E03-MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A01-FILTER_MUTECT2_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${SM_TAG} \
					${REF_GENOME} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		###################################################
		# apply masks to mutect2 mito filtered vcf output #
		###################################################

			MASK_MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A01-A01-MASK_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-MASK_MUTECT2_MT.log \
				-hold_jid E03-A01-FILTER_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A01-A01-MASK_MUTECT2_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${MT_MASK} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		#############################################
		# run haplogrep2 on mutect2 mito vcf output #
		#############################################

			HAPLOGREP2_MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A01-A01-A01-HAPLOGREP2_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-HAPLOGREP2_MUTECT2_MT.log \
				-hold_jid E03-A01-A01-MASK_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A01-A01-A01-HAPLOGREP2_MUTECT2_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${REF_GENOME} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		##########################################
		# run annovar on final mutect2 based vcf #
		##########################################

			RUN_ANNOVAR_MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A01-A01-A02-RUN_ANNOVAR_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_ANNOVAR_MUTECT2_MT.log \
				-hold_jid E03-A01-A01-MASK_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A01-A01-A02-RUN_ANNOVAR_MUTECT2_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${ANNOVAR_MT_DB_DIR} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		##########################################
		# run annovar on final mutect2 based vcf #
		##########################################

			FIX_ANNOVAR_MUTECT2_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A01-A01-A02-A01-FIX_ANNOVAR_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FIX_ANNOVAR_MUTECT2_MT.log \
				-hold_jid E03-A01-A01-A02-RUN_ANNOVAR_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A01-A01-A02-A01-FIX_ANNOVAR_MUTECT2_MT.sh \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		##################################
		# CONVERT MUTECT2 MT BAM TO CRAM #
		##################################

			MUTECT2_MT_BAM_TO_CRAM ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E03-A02-MUTECT2_MT_BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-MUTECT2_MT_BAM_TO_CRAM.log \
				-hold_jid E03-MUTECT2_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E03-A02-MUTECT2_MT_BAM_TO_CRAM.sh \
					${ALIGNMENT_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${REF_GENOME} \
					${THREADS} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

	###############################################################
	### EKLIPSE WORKFLOW TO DETECT LARGE DELETIONS IN MT GENOME ###
	###############################################################

		############################################
		# SUBSET BAM FILE TO CONTAIN ONLY MT READS #
		############################################

			SUBSET_BAM_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E04-MAKE_BAM_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-MAKE_BAM_MT.log \
				-hold_jid D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E04-MAKE_MT_BAM.sh \
					${MITO_EKLIPSE_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${THREADS} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		###############
		# RUN EKLIPSE #
		###############

			RUN_EKLIPSE ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E04-A01-RUN_EKLIPSE_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_EKLIPSE.log \
				-hold_jid E04-MAKE_BAM_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E04-A01-RUN_EKLIPSE.sh \
					${MITO_EKLIPSE_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${MT_GENBANK} \
					${THREADS} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		#####################################
		# ADD LEGEND TO EKLIPSE CIRCOS PLOT #
		#####################################

			FORMAT_EKLIPSE_CIRCOS ()
			{
				echo \
				qsub \
					${IMGMAGICK_QSUB_ARGS} \
				-N E04-A01-A01-FORMAT_EKLIPSE_CIRCOS_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FORMAT_EKLIPSE_CIRCOS.log \
				-hold_jid E04-A01-RUN_EKLIPSE_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E04-A01-A01-FORMAT_EKLIPSE_CIRCOS.sh \
					${MITO_MAGICK_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${EKLIPSE_FORMAT_CIRCOS_PLOT_R_SCRIPT} \
					${EKLIPSE_CIRCOS_LEGEND} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

	##################################################
	### COVERAGE STATISTICS AND PLOT FOR MT GENOME ###
	##################################################

		####################################################
		# RUN COLLECTHSMETRICS ON MT ONLY BAM FILE #########
		# USES GATK IMPLEMENTATION INSTEAD OF PICARD TOOLS #
		####################################################

			COLLECTHSMETRICS_MT ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E04-A02-COLLECTHSMETRICS_MT_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-COLLECTHSMETRICS_MT.log \
				-hold_jid E04-MAKE_BAM_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E04-A02-COLLECTHSMETRICS_MT.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${REF_GENOME} \
					${MT_PICARD_INTERVAL_LIST} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

		###############################################################
		# RUN ALEX'S R SCRIPT TO GENERATE COVERAGE PLOT FOR MT GENOME #
		###############################################################

			PLOT_MT_COVERAGE ()
			{
				echo \
				qsub \
					${QSUB_ARGS} \
				-N E04-A02-A01-PLOT_MT_COVERAGE_${SGE_SM_TAG}_${PROJECT} \
					-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-PLOT_MT_COVERAGE.log \
				-hold_jid E04-A02-COLLECTHSMETRICS_MT_${SGE_SM_TAG}_${PROJECT} \
				${SCRIPT_DIR}/E04-A02-A01-PLOT_MT_COVERAGE.sh \
					${MITO_MUTECT2_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${MT_COVERAGE_R_SCRIPT} \
					${SAMPLE_SHEET} \
					${SUBMIT_STAMP}
			}

############################################
# RUN STEPS FOR MITOCHONDRIAL DNA ANALYSIS #
############################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
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
		RUN_ANNOVAR_MUTECT2_MT
		echo sleep 0.1s
		FIX_ANNOVAR_MUTECT2_MT
		echo sleep 0.1s
		MUTECT2_MT_BAM_TO_CRAM
		echo sleep 0.1s
		# run eklipse workflow
		SUBSET_BAM_MT
		echo sleep 0.1s
		RUN_EKLIPSE
		echo sleep 0.1s
		FORMAT_EKLIPSE_CIRCOS
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
				${QSUB_ARGS} \
			-N E05-RUN_EXOME_DEPTH_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_EXOME_DEPTH.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E05-RUN_EXOME_DEPTH.sh \
				${CNV_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${EXOME_DEPTH_R_SCRIPT} \
				${REF_PANEL_COUNTS} \
				${CODING_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	################################################################
	# calculate the percent of cnv call length for each chromosome #
	################################################################

		CALCULATE_PCT_CNV_COVERAGE ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E05-A01-PCT_CNV_COVERAGE_PER_CHR_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-PCT_CNV_COVERAGE_PER_CHR.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E05-RUN_EXOME_DEPTH_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E05-A01-PCT_CNV_COVERAGE_PER_CHR.sh \
				${CNV_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED}
		}

	########################################
	# run annotSV on the exomeDepth output #
	########################################

		RUN_ANNOTSV ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E05-A02-RUN_ANNOTSV_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_ANNOTSV.log \
			-hold_jid E05-RUN_EXOME_DEPTH_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E05-A02-RUN_ANNOTSV.sh \
				${CNV_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	########################################################################################
	# reformat the header in the annotSV output and filter to zoom gene list if applicable #
	########################################################################################

		RUN_FORMAT_AND_ZOOM_ANNOTSV ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E05-A02-A01-RUN_FORMAT_AND_ZOOM_ANNOTSV_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_FORMAT_AND_ZOOM_ANNOTSV.log \
			-hold_jid E05-A02-RUN_ANNOTSV_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E05-A02-A01-RUN_FORMAT_AND_ZOOM_ANNOTSV.sh \
				${CNV_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${FORMAT_AND_ZOOM_ANNOTSV_R_SCRIPT} \
				${ZOOM_LIST} \
				${ZOOM_NAME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

##############################
# RUN STEPS FOR CNV WORKFLOW #
##############################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
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

####################################################################################
##### BAM FILE METRICS #############################################################
####################################################################################
# NOTE: SOME PROGRAMS CAN ONLY BE RAN ON THE BAM FILE AND NOT ON THE CRAM FILE #####
####################################################################################

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR CODING BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container ##################################################
	# This with all RefSeq Select CDS exons + missing OMIM, etc. #################
	##############################################################################

		DOC_CODING ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-DOC_CODING_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-DOC_CODING.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-DOC_CODING_PADDED.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${CODING_BED} \
				${PADDING_LENGTH} \
				${GENE_LIST} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#########################################################
	# DO AN ANEUPLOIDY CHECK ON CODING BED FILE DOC OUTPUT  #
	#########################################################

		ANEUPLOIDY_CHECK ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A01-CHROM_DEPTH_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-ANEUPLOIDY_CHECK.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E06-DOC_CODING_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A01-CHROM_DEPTH.sh \
					${ALIGNMENT_CONTAINER} \
					${CORE_PATH} \
					${PROJECT} \
					${FAMILY} \
					${SM_TAG} \
					${CODING_BED} \
					${PADDING_LENGTH}
		}

	########################################################################################
	# FORMATTING PER BASE COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	########################################################################################

		ANNOTATE_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A02-ANNOTATE_PER_BASE_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-ANNOTATE_PER_BASE.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E06-DOC_CODING_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A02-ANNOTATE_PER_BASE.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED} \
				${PADDING_LENGTH} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	##########################################################################
	# FILTER PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x #
	##########################################################################

		FILTER_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A02-A01-FILTER_ANNOTATED_PER_BASE_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_ANNOTATED_PER_BASE.log \
			-hold_jid E06-A02-ANNOTATE_PER_BASE_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A02-A01-FILTER_ANNOTATED_PER_BASE.sh \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED} \
				${PADDING_LENGTH}
		}

	######################################################
	# BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION #
	######################################################

		BGZIP_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A02-A01-A01-BGZIP_ANNOTATED_PER_BASE_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-BGZIP_ANNOTATED_PER_BASE.log \
			-hold_jid E06-A02-A01-FILTER_ANNOTATED_PER_BASE_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A02-A01-A01-BGZIP_ANNOTATED_PER_BASE.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED} \
				${PADDING_LENGTH} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	###################################################################################################
	# FORMATTING PER CODING INTERVAL COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	###################################################################################################

		ANNOTATE_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A03-ANNOTATE_PER_INTERVAL_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-ANNOTATE_PER_INTERVAL.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E06-DOC_CODING_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A03-ANNOTATE_PER_INTERVAL.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED} \
				${PADDING_LENGTH}
		}

	##################################################################################################
	# FILTER ANNOTATED PER CODING INTERVAL COVERAGE TO INTERVALS WHERE LESS 100% OF BASES ARE AT 30X #
	##################################################################################################

		FILTER_ANNOTATED_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E06-A03-A01-FILTER_ANNOTATED_PER_INTERVAL_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_ANNOTATED_PER_INTERVAL.log \
			-hold_jid E06-A03-ANNOTATE_PER_INTERVAL_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E06-A03-A01-FILTER_ANNOTATED_PER_INTERVAL.sh \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${CODING_BED} \
				${PADDING_LENGTH}
		}

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR TARGET BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container #############################################################################
	# Generally this with all RefSeq Select CDS exons + missing OMIM unless it becomes targeted, e.g a zoom #
	#########################################################################################################

		DOC_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E07-DOC_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-DOC_TARGET.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E07-DOC_TARGET_PADDED.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${TARGET_BED} \
				${PADDING_LENGTH} \
				${GENE_LIST} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
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
				${QSUB_ARGS} \
			-N E08-SELECT_VERIFYBAMID_VCF_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-SELECT_VERIFYBAMID_VCF.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E08-SELECT_VERIFYBAMID_VCF.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${SM_TAG} \
				${REF_GENOME} \
				${VERIFY_VCF} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	###################
	# RUN VERIFYBAMID #
	###################

		RUN_VERIFYBAMID ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E08-A01-RUN_VERIFYBAMID_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-VERIFYBAMID.log \
			-hold_jid E08-SELECT_VERIFYBAMID_VCF_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E08-A01-VERIFYBAMID.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#################################################################################################
	# CREATE VCF PER CHROMOSOME AND RUN VERIFYBAMID ON THEM ######################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	# SCRIPT READS BAIT BED FILE, GRABS THE CHROMOSOMES AND RUNS A FOR LOOP FOR BOTH THINGS #########
	#################################################################################################

		VERIFYBAMID_PER_AUTOSOME ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N E09-SELECT_VERIFYBAMID_PER_AUTOSOME_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-SELECT_VERIFYBAMID_PER_AUTOSOME.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},D01-APPLY_BQSR_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E09-VERIFYBAMID_PER_AUTO.sh \
				${ALIGNMENT_CONTAINER} \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${SM_TAG} \
				${REF_GENOME} \
				${VERIFY_VCF} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
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
				${QSUB_ARGS} \
			-N E09-A01-CAT_VERIFYBAMID_AUTOSOME_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-CAT_VERIFYBAMID_AUTOSOME.log \
			-hold_jid E09-SELECT_VERIFYBAMID_PER_AUTOSOME_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/E09-A01-CAT_VERIFYBAMID_AUTO.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${BAIT_BED}
		}

##########################################
# RUN STEPS FOR BAM FILE RELATED METRICS #
##########################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_SAMPLE_ARRAY
		DOC_CODING
		echo sleep 0.1s
		DOC_TARGET
		echo sleep 0.1s
		SELECT_VERIFYBAMID_VCF
		echo sleep 0.1s
		RUN_VERIFYBAMID
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

##################################################################
##### JOINT CALLING SAMPLES IN A FAMILY WITH SET OF CONTROLS #####
##################################################################

	###############################################################################################
	# create a hold_id variable for haplotype caller gvcf gather step for all samples in a family #
	###############################################################################################

		BUILD_HOLD_ID_PATH_GENOTYPE_GVCF ()
		{
			for PROJECT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					{print $1}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq )
			do
				HOLD_ID_PATH="-hold_jid "

				for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
						$20=="'${FAMILY_ONLY}'" \
						{print $8}' \
					~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
						| sed 's/@/_/g' \
						| sort \
						| uniq);
				do
					HOLD_ID_PATH="${HOLD_ID_PATH}E02-A01-HAPLOTYPE_CALLER_GVCF_GATHER_${SAMPLE}_${PROJECT},"
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
				${QSUB_ARGS} \
			-N F01-GENOTYPE_GVCF_SCATTER_${FAMILY}_${PROJECT}_chr${CHROMOSOME} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.GENOTYPE_GVCF_chr${CHROMOSOME}.log \
			${HOLD_ID_PATH} \
			${SCRIPT_DIR}/F01-GENOTYPE_GVCF_SCATTER.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${DBSNP} \
				${CHROMOSOME} \
				${CONTROL_REPO} \
				${CONTROL_DATA_SET_FILE} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	########################################
	# scatter genotype gvcfs by chromosome #
	########################################

		SCATTER_GENOTYPE_GVCF_PER_CHROMOSOME ()
		{
			for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/chr//g' \
				| egrep "^[0-9]|^X|^Y" \
				| cut -f 1 \
				| sort -V \
				| uniq \
				| singularity exec ${ALIGNMENT_CONTAINER} datamash \
					collapse 1 \
				| sed 's/,/ /g');
			do
				GENOTYPE_GVCF
				echo sleep 0.1s
			done
		}

#################################################################################
# RUN STEPS TO DO JOINT CALLING PER FAMILY PER SET OF INTERVALS IN A CHROMOSOME #
#################################################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_FAMILY_ARRAY
		BUILD_HOLD_ID_PATH_GENOTYPE_GVCF
		SCATTER_GENOTYPE_GVCF_PER_CHROMOSOME
		echo sleep 0.1s
	done

#############################
##### CRAM FILE METRICS #####
#############################

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
				${QSUB_ARGS} \
			-N F02-COLLECT_MULTIPLE_METRICS_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-COLLECT_MULTIPLE_METRICS.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/F02-COLLECT_MULTIPLE_METRICS.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${DBSNP} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
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
				${QSUB_ARGS} \
			-N F03-COLLECT_HS_METRICS_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-COLLECT_HS_METRICS.log \
			-hold_jid A02-FIX_BED_FILES_${SGE_SM_TAG}_${PROJECT},E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/F03-COLLECT_HS_METRICS.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${BAIT_BED} \
				${TITV_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

###########################################
# RUN STEPS FOR CRAM FILE RELATED METRICS #
###########################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_SAMPLE_ARRAY
		COLLECT_MULTIPLE_METRICS
		echo sleep 0.1s
		COLLECT_HS_METRICS
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
			for PROJECT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					{print $1}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq )
			do
				HOLD_ID_PATH="-hold_jid "

				for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
					| sed -r 's/[[:space:]]+/\t/g' \
					| sed 's/chr//g' \
					| egrep "^[0-9]|^X|^Y" \
					| cut -f 1 \
					| sort -V \
					| uniq \
					| singularity exec ${ALIGNMENT_CONTAINER} datamash \
						collapse 1 \
					| sed 's/,/ /g');
				do
					HOLD_ID_PATH="${HOLD_ID_PATH}F01-GENOTYPE_GVCF_SCATTER_${FAMILY}_${PROJECT}_chr${CHROMOSOME},"
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
				${QSUB_ARGS} \
			-N G01-GENOTYPE_GVCF_GATHER_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.GENOTYPE_GVCF_GATHER.log \
			${HOLD_ID_PATH}A03-FIX_BED_FILES_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/G01-GENOTYPE_GVCF_GATHER.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

#####################################################
# RUN STEP TO GATHER PER CHROMOSOME PER FAMILY VCFS #
#####################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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
# BUT I LIKE TO KEEP IT HERE ###########################
########################################################

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		RUN_VQSR_SNP ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N H01-RUN_VQSR_SNP_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.RUN_VQSR_SNP.log \
			-hold_jid G01-GENOTYPE_GVCF_GATHER_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/H01-RUN_VARIANT_RECALIBRATOR_SNP.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${DBSNP} \
				${HAPMAP} \
				${OMNI_1KG} \
				${HI_CONF_1KG_PHASE1_SNP} \
				${SEND_TO} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	################################################
	# Run Variant Recalibrator for the INDEL model #
	################################################

		RUN_VQSR_INDEL ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N H02-RUN_VQSR_INDEL_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.RUN_VQSR_INDEL.log \
			-hold_jid G01-GENOTYPE_GVCF_GATHER_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/H02-RUN_VARIANT_RECALIBRATOR_INDEL.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${MILLS_1KG_GOLD_INDEL} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		APPLY_VQSR_SNP ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N I01-APPLY_VQSR_SNP_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.APPLY_VQSR_SNP.log \
			-hold_jid H01-RUN_VQSR_SNP_${FAMILY}_${PROJECT},H02-RUN_VQSR_INDEL_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/I01-APPLY_VARIANT_RECALIBRATION_SNP.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	##############################################
	# Run Variant Recalibrator for the SNP model #
	##############################################

		APPLY_VQSR_INDEL ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N J01-APPLY_VQSR_INDEL_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.APPLY_VQSR_INDEL.log \
			-hold_jid I01-APPLY_VQSR_SNP_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/J01-APPLY_VARIANT_RECALIBRATION_INDEL.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

########################
# RUN STEPS TO DO VQSR #
########################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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
			${QSUB_ARGS} \
		-N K01-VARIANT_ANNOTATOR_${FAMILY}_${PROJECT}_${CHROMOSOME} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.VARIANT_ANNOTATOR_${CHROMOSOME}.log \
		-hold_jid J01-APPLY_VQSR_INDEL_${FAMILY}_${PROJECT} \
		${SCRIPT_DIR}/K01-VARIANT_ANNOTATOR_SCATTER.sh \
			${GATK_3_7_0_CONTAINER} \
			${CORE_PATH} \
			${PED_FILE} \
			${PROJECT} \
			${FAMILY} \
			${REF_GENOME} \
			${CHROMOSOME} \
			${PHASE3_1KG_AUTOSOMES} \
			${THREADS} \
			${SAMPLE_SHEET} \
			${SUBMIT_STAMP}
	}

#####################################
# RUN STEPS TO DO VARIANT ANNOTATOR #
#####################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_FAMILY_ARRAY

		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
			| sed -r 's/[[:space:]]+/\t/g' \
			| sed 's/chr//g' \
			| egrep "^[0-9]|^X|^Y" \
			| cut -f 1 \
			| sort -V \
			| uniq \
			| singularity exec ${ALIGNMENT_CONTAINER} datamash \
				collapse 1 \
			| sed 's/,/ /g');
		do
			CALL_VARIANT_ANNOTATOR
			echo sleep 0.1s
		done
	done

##############################################################################################
##### GATHER UP THE PER FAMILY PER CHROMOSOME ANNOTATED VCF FILES INTO A SINGLE VCF FILE #####
##### RUN PCA/RELATEDNESS WORKFLOW ###########################################################
##### RUN ROH ANALYSIS #######################################################################
##############################################################################################

	######################################################
	# generate hold id from scatter of variant annotator #
	######################################################

		BUILD_HOLD_ID_PATH_ADD_MORE_ANNOTATION ()
		{
			for PROJECT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					{print $1}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq )
			do
				HOLD_ID_PATH="-hold_jid "

				for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
					| sed -r 's/[[:space:]]+/\t/g' \
					| sed 's/chr//g' \
					| egrep "^[0-9]|^X|^Y" \
					| cut -f 1 \
					| sort -V \
					| uniq \
					| singularity exec ${ALIGNMENT_CONTAINER} datamash \
						collapse 1 \
					| sed 's/,/ /g');
				do
					HOLD_ID_PATH="${HOLD_ID_PATH}K01-VARIANT_ANNOTATOR_${FAMILY}_${PROJECT}_${CHROMOSOME},"
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
				${QSUB_ARGS} \
			-N L01-VARIANT_ANNOTATOR_GATHER_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.MORE_VARIANT_ANNOTATOR_GATHER.log \
			${HOLD_ID_PATH} \
			${SCRIPT_DIR}/L01-VARIANT_ANNOTATOR_GATHER.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#####################################################################
	# FILTER TO JUST PASSING BIALLELIC SNV SITES ON THE CODING BED FILE #
	# TEMPORARY FILE USED FOR PCA AND RELATEDNESS #######################
	#####################################################################

		CALL_PASS_BIALLELIC_SNV_COHORT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N M01-FILTER_COHORT_SNV_PASS_BIALLELIC_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_COHORT_SNV_PASS_BIALLELIC.log \
			-hold_jid L01-VARIANT_ANNOTATOR_GATHER_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/M01-FILTER_COHORT_SNV_PASS_BIALLELIC.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${CODING_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	################################
	# RUN PCA AND KINSHIP WORKFLOW #
	# USES KING AND PLINK ##########
	################################

		CALL_PCA_RELATEDNESS ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N M01-A01-PCA_RELATEDNESS_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.PCA_RELATEDNESS.log \
			-hold_jid M01-FILTER_COHORT_SNV_PASS_BIALLELIC_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/M01-A01-PCA_RELATEDNESS.sh \
				${GATK_3_7_0_CONTAINER} \
				${PCA_RELATEDNESS_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${PED_FILE} \
				${CONTROL_PED_FILE} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	###########################################################
	# FILTER OUT REPEATMASKED REGIONS TO PERFORM ROH ANALYSIS #
	###########################################################

		FILTER_REPEATMASK ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N M01-A02-FILTER_REPEATMASK_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_REPEATMASK.log \
			-hold_jid M01-FILTER_COHORT_SNV_PASS_BIALLELIC_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/M01-A02-FILTER_REPEATMASK.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${UCSC_REPEATMASK} \
				${MDUST_REPEATMASK} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	####################
	# RUN ROH ANALYSIS #
	####################

		BCFTOOLS_ROH ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N M01-A02-A01-BCFTOOLS_ROH_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.BCFTOOLS_ROH.log \
			-hold_jid M01-A02-FILTER_REPEATMASK_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/M01-A02-A01-BCFTOOLS_ROH.sh \
				${MITO_MUTECT2_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

##################################################################
# RUN STEPS TO DO VARIANT ANNOTATOR GATHER, PCA/RELATEDNESS, ROH #
##################################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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
		FILTER_REPEATMASK
		echo sleep 0.1s
		BCFTOOLS_ROH
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
			${QSUB_ARGS} \
		-N L02-FILTER_TO_FAMILY_ALL_SITES_${FAMILY}_${PROJECT}_${CHROMOSOME} \
			-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_TO_FAMILY_ALL_SITES_${CHROMOSOME}.log \
		-hold_jid K01-VARIANT_ANNOTATOR_${FAMILY}_${PROJECT}_${CHROMOSOME} \
		${SCRIPT_DIR}/L02-FILTER_TO_FAMILY_ALL_SITES_CHR.sh \
			${ALIGNMENT_CONTAINER} \
			${CORE_PATH} \
			${PROJECT} \
			${FAMILY} \
			${REF_GENOME} \
			${CHROMOSOME} \
			${SAMPLE_SHEET} \
			${SUBMIT_STAMP}
	}

####################################################
# RUN STEPS TO FILTER ALL SITES VCF TO FAMILY ONLY #
####################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
	do
		CREATE_FAMILY_ARRAY

		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
			| sed -r 's/[[:space:]]+/\t/g' \
			| sed 's/chr//g' \
			| egrep "^[0-9]|^X|^Y" \
			| cut -f 1 \
			| sort -V \
			| uniq \
			| singularity exec ${ALIGNMENT_CONTAINER} datamash \
				collapse 1 \
			| sed 's/,/ /g');
		do
			CALL_FILTER_TO_FAMILY_ALL_SITES
			echo sleep 0.1s
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
			for PROJECT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					{print $1}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq );
			do
				HOLD_ID_PATH="-hold_jid "

				for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' ${BAIT_BED} \
					| sed -r 's/[[:space:]]+/\t/g' \
					| sed 's/chr//g' \
					| egrep "^[0-9]|^X|^Y" \
					| cut -f 1 \
					| sort -V \
					| uniq \
					| singularity exec ${ALIGNMENT_CONTAINER} datamash \
						collapse 1 \
					| sed 's/,/ /g');
				do
					HOLD_ID_PATH="${HOLD_ID_PATH}L02-FILTER_TO_FAMILY_ALL_SITES_${FAMILY}_${PROJECT}_${CHROMOSOME},"
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
				${QSUB_ARGS} \
			-N M02-FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_TO_FAMILY_ALL_SITES_GATHER.log \
			${HOLD_ID_PATH} \
			${SCRIPT_DIR}/M02-FILTER_TO_FAMILY_ALL_SITES_GATHER.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${BAIT_BED} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	############################################################
	# filter family only all sites vcf to coding plus user pad #
	# the output for this might get moved to temp ##############
	############################################################

		CALL_FILTER_FAMILY_TO_CODING_PLUS_PAD ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N N01-FILTER_FAMILY_CODING_PLUS_PAD_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_FAMILY_CODING_PLUS_PAD.log \
			-hold_jid M02-FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/N01-FILTER_FAMILY_CODING_PLUS_PAD.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${CODING_BED} \
				${PADDING_LENGTH} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	############################################################
	# filter family only all sites vcf to target plus user pad #
	############################################################

		CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N N02-FILTER_FAMILY_TARGET_PLUS_PAD_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_FAMILY_TARGET_PLUS_PAD.log \
			-hold_jid M02-FILTER_TO_FAMILY_ALL_SITES_GATHER_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/N02-FILTER_FAMILY_TARGET_PLUS_PAD.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${TARGET_BED} \
				${PADDING_LENGTH} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	############################################################
	# filter family only all sites vcf to target plus user pad #
	############################################################

		CALL_FILTER_FAMILY_TO_TARGET_PLUS_PAD_VARIANTS ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N N02-A01-FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS_${FAMILY}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/LOGS/${FAMILY}_${PROJECT}.FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS.log \
			-hold_jid N02-FILTER_FAMILY_TARGET_PLUS_PAD_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/N02-A01-FILTER_TO_FAMILY_TARGET_PLUS_PAD_VARIANTS.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

########################################################################
# RUN STEPS TO GATHER UP PER CHROMOSOME FAMILY ONLY ALL SITES VCF FILE #
########################################################################

	for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $20}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
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

#####################################
##### SUBSETTING TO SAMPLE VCFS #####
#####################################

	#####################################################################################
	# subset sample all sites to from family coding/bait bed file plus user defined pad #
	#####################################################################################

		EXTRACT_SAMPLE_ALL_SITES ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N P01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_ALL_SITES.log \
			-hold_jid N01-FILTER_FAMILY_CODING_PLUS_PAD_${FAMILY}_${PROJECT} \
			${SCRIPT_DIR}/P01-FILTER_TO_SAMPLE_ALL_SITES.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	##################################################################################
	# subset sample variant sites to from coding/bait bed file plus user defined pad #
	##################################################################################

		EXTRACT_SAMPLE_VARIANTS ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N Q01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_VARIANTS.log \
			-hold_jid P01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/Q01-FILTER_TO_SAMPLE_VARIANTS.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#################################################################################################
	# generate vcf metrics for sample variant sites from coding/bait bed file plus user defined pad #
	#################################################################################################

		VCF_METRICS_BAIT ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N Q01-A01-VCF_METRICS_BAIT_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-VCF_METRICS_BAIT.log \
			-hold_jid Q01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/Q01-A01-VCF_METRICS_BAIT.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_DICT} \
				${DBSNP} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#####################################################################
	# generate vcf metrics for sample variant sites from ti/tv bed file #
	#####################################################################

		VCF_METRICS_TITV ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N Q01-A02-VCF_METRICS_TITV_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-VCF_METRICS_TITV.log \
			-hold_jid Q01-FILTER_TO_SAMPLE_VARIANTS_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/Q01-A02-VCF_METRICS_TITV.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_DICT} \
				${TITV_BED} \
				${DBSNP_129} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#########################################################################
	# subset sample to all sites from target bed file plus user defined pad #
	#########################################################################

		EXTRACT_SAMPLE_ALL_SITES_ON_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N Q02-FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_ALL_SITES_TARGET.log \
			-hold_jid P01-FILTER_TO_SAMPLE_ALL_SITES_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/Q02-FILTER_TO_SAMPLE_ALL_SITES_TARGET.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${TARGET_BED} \
				${PADDING_LENGTH} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	########################################################################
	# subset sample variant sites to target bed file plus user defined pad #
	########################################################################

		EXTRACT_SAMPLE_VARIANTS_ON_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N R01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-FILTER_TO_SAMPLE_VARIANTS_TARGET.log \
			-hold_jid Q02-FILTER_TO_SAMPLE_ALL_SITES_TARGET_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/R01-FILTER_TO_SAMPLE_VARIANTS_TARGET.sh \
				${GATK_3_7_0_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#####################################################################################
	# decompose on target plus user defined pad sample variant sites to target bed file #
	#####################################################################################

		DECOMPOSE_SAMPLE_VARIANTS_ON_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N R01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-DECOMPOSE_SAMPLE_VARIANTS_TARGET.log \
			-hold_jid R01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/R01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET.sh \
				${VT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_GENOME} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	#############################################################################
	# run annovar on decomposed on target plus user defined pad sample vcf file #
	#############################################################################

		RUN_ANNOVAR_ON_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N R01-A01-A01-RUN_ANNOVAR_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-RUN_ANNOVAR_TARGET.log \
			-hold_jid R01-A01-DECOMPOSE_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/R01-A01-A01-RUN_ANNOVAR_TARGET.sh \
				${ANNOVAR_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${ANNOVAR_DATABASE_FILE} \
				${ANNOVAR_REF_BUILD} \
				${ANNOVAR_INFO_FIELD_KEYS} \
				${ANNOVAR_HEADER_MAPPINGS} \
				${ANNOVAR_VCF_COLUMNS} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

	############################################################################################
	# generate vcf metrics for sample variant sites from target bed file plus user defined pad #
	############################################################################################

		VCF_METRICS_TARGET ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N R01-A02-VCF_METRICS_TARGET_${SGE_SM_TAG}_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-VCF_METRICS_TARGET.log \
			-hold_jid R01-FILTER_TO_SAMPLE_VARIANTS_TARGET_${SGE_SM_TAG}_${PROJECT} \
			${SCRIPT_DIR}/R01-A02-VCF_METRICS_TARGET.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${FAMILY} \
				${SM_TAG} \
				${REF_DICT} \
				${DBSNP} \
				${THREADS} \
				${SAMPLE_SHEET} \
				${SUBMIT_STAMP}
		}

######################################
### QC REPORT PREP FOR EACH SAMPLE ###
######################################

QC_REPORT_PREP ()
{
echo \
qsub \
	${QSUB_ARGS} \
-N X01-QC_REPORT_PREP_${SGE_SM_TAG}_${PROJECT} \
	-o ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/LOGS/${SM_TAG}-QC_REPORT_PREP.log \
-hold_jid \
R01-A02-VCF_METRICS_TARGET_${SGE_SM_TAG}_${PROJECT},\
R01-A01-A01-RUN_ANNOVAR_TARGET_${SGE_SM_TAG}_${PROJECT},\
Q01-A02-VCF_METRICS_TITV_${SGE_SM_TAG}_${PROJECT},\
Q01-A01-VCF_METRICS_BAIT_${SGE_SM_TAG}_${PROJECT},\
E03-A02-MUTECT2_MT_BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT},\
E04-A02-A01-PLOT_MT_COVERAGE_${SGE_SM_TAG}_${PROJECT},\
E04-A01-A01-FORMAT_EKLIPSE_CIRCOS_${SGE_SM_TAG}_${PROJECT},\
E03-A01-A01-A02-A01-FIX_ANNOVAR_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT},\
E03-A01-A01-A01-HAPLOGREP2_MUTECT2_MT_${SGE_SM_TAG}_${PROJECT},\
E09-A01-CAT_VERIFYBAMID_AUTOSOME_${SGE_SM_TAG}_${PROJECT},\
E06-A03-A01-FILTER_ANNOTATED_PER_INTERVAL_${SGE_SM_TAG}_${PROJECT},\
E06-A02-A01-A01-BGZIP_ANNOTATED_PER_BASE_${SGE_SM_TAG}_${PROJECT},\
E06-A02-A01-FILTER_ANNOTATED_PER_BASE_${SGE_SM_TAG}_${PROJECT},\
E06-A01-CHROM_DEPTH_${SGE_SM_TAG}_${PROJECT},\
E08-A01-RUN_VERIFYBAMID_${SGE_SM_TAG}_${PROJECT},\
E07-DOC_TARGET_${SGE_SM_TAG}_${PROJECT},\
F03-COLLECT_HS_METRICS_${SGE_SM_TAG}_${PROJECT},\
F02-COLLECT_MULTIPLE_METRICS_${SGE_SM_TAG}_${PROJECT},\
E05-A01-PCT_CNV_COVERAGE_PER_CHR_${SGE_SM_TAG}_${PROJECT},\
E05-A02-A01-RUN_FORMAT_AND_ZOOM_ANNOTSV_${SGE_SM_TAG}_${PROJECT},\
E01-BAM_TO_CRAM_${SGE_SM_TAG}_${PROJECT} \
${SCRIPT_DIR}/X01-QC_REPORT_PREP.sh \
	${ALIGNMENT_CONTAINER} \
	${CORE_PATH} \
	${PROJECT} \
	${FAMILY} \
	${SM_TAG} \
	${FATHER} \
	${MOTHER} \
	${GENDER} \
	${PHENOTYPE}
}

##################################################
# RUN VCF SAMPLE SUBSTEP STEP AND QC REPORT PREP #
##################################################

	for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $8}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq);
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

#############################
##### END PROJECT TASKS #####
#############################

	############################################################
	# build hold id for qc report prep per sample, per project #
	############################################################

		BUILD_HOLD_ID_PATH_PROJECT_WRAP_UP_SAMPLE ()
		{
			HOLD_ID_PATH_QC_REPORT_PREP="-hold_jid "

			for SAMPLE in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					$1=="'${PROJECT}'" \
					{print $8}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq);
			do
				CREATE_SAMPLE_ARRAY

				HOLD_ID_PATH_QC_REPORT_PREP="${HOLD_ID_PATH_QC_REPORT_PREP}X01-QC_REPORT_PREP_${SGE_SM_TAG}_${PROJECT},"

				HOLD_ID_PATH_QC_REPORT_PREP=`echo ${HOLD_ID_PATH_QC_REPORT_PREP} | sed 's/@/_/g'`
			done
		}

	###################################################
	# add hold id for PCA/RELATEDNESS FOR EACH FAMILY #
	###################################################

		BUILD_HOLD_ID_PATH_PROJECT_WRAP_UP_FAMILY ()
		{
			HOLD_ID_PATH_PCA=""

			for FAMILY_ONLY in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
					$1=="'${PROJECT}'" \
					{print $20}' \
				~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
					| sort \
					| uniq);
			do
				CREATE_FAMILY_ARRAY

				HOLD_ID_PATH_PCA="${HOLD_ID_PATH_PCA}M01-A01-PCA_RELATEDNESS_${FAMILY}_${PROJECT},M01-A02-A01-BCFTOOLS_ROH_${FAMILY}_${PROJECT},"
			done
		}

	#########################################################################
	# run end project functions (qc report, file clean-up) for each project #
	#########################################################################

		PROJECT_WRAP_UP ()
		{
			echo \
			qsub \
				${QSUB_ARGS} \
			-N X01-X01_END_PROJECT_TASKS_${PROJECT} \
				-o ${CORE_PATH}/${PROJECT}/LOGS/${PROJECT}-END_PROJECT_TASKS.log \
			${HOLD_ID_PATH_QC_REPORT_PREP}${HOLD_ID_PATH_PCA} \
			${SCRIPT_DIR}/X01-X01-END_PROJECT_TASKS.sh \
				${ALIGNMENT_CONTAINER} \
				${CORE_PATH} \
				${PROJECT} \
				${SCRIPT_DIR} \
				$SUBMITTER_ID \
				${SAMPLE_SHEET} \
				${PED_FILE} \
				${SUBMIT_STAMP} \
				${SEND_TO} \
				${THREADS}
		}

##################
# RUN FINAL LOOP #
##################

	for PROJECT in $(awk 'BEGIN {FS="\t"; OFS="\t"} \
			{print $1}' \
		~/CGC_PIPELINE_TEMP/${MANIFEST_PREFIX}.${PED_PREFIX}.join.txt \
			| sort \
			| uniq );
	do
		BUILD_HOLD_ID_PATH_PROJECT_WRAP_UP_SAMPLE
		BUILD_HOLD_ID_PATH_PROJECT_WRAP_UP_FAMILY
		PROJECT_WRAP_UP
	done

#############################################################
##### MESSAGE THAT SAMPLE SHEET HAS FINISHED SUBMITTING #####
#############################################################

	printf "echo\n"

	printf "echo ${SAMPLE_SHEET} has finished submitting at `date`\n"

######################################
##### EMAIL WHEN DONE SUBMITTING #####
######################################

	printf "${SAMPLE_SHEET}\nhas finished submitting at\n`date`\nby `whoami`" \
		| mail -s "${PERSON_NAME} has submitted SUBMITTER_JHG-DDL_EXOME_PIPELINE.sh" \
			${SEND_TO}
