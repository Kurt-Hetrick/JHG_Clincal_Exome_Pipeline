#!/bin/bash

JAVA_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jdk1.8.0_73/bin"
GATK_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/GenomeAnalysisTK-3.7"
REFERENCE_GENOME="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/human_g1k_v37_decoy.fasta"

GVCF_LIST=$1
OUTPUT=$2

$JAVA_DIR/java -jar \
$GATK_DIR/GenomeAnalysisTK.jar \
-T CombineGVCFs \
--reference_sequence $REFERENCE_GENOME \
--annotation AS_BaseQualityRankSumTest \
--annotation AS_FisherStrand \
--annotation AS_InbreedingCoeff \
--annotation AS_MappingQualityRankSumTest \
--annotation AS_RMSMappingQuality \
--annotation AS_ReadPosRankSumTest \
--annotation AS_StrandOddsRatio \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--variant $GVCF_LIST \
--out $OUTPUT
