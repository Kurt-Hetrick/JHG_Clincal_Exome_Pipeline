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

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3
VCFTOOLS_DIR=$4
PLINK2_DIR=$5
KING_DIR=$6

PROJECT=$7
FAMILY=$8
REF_GENOME=$9
PED_FILE=${10}
CONTROL_PED_FILE=${11}

# Format the control ped file in case molly or somebody ever changes it and screws up the format

awk 1 $CONTROL_PED_FILE \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
>| $CORE_PATH/$PROJECT/TEMP/"CONTROL_PED_FILE_FOR_"$FAMILY".ped"

# Concatenate the control ped file with the ped information for the family

awk 1 $PED_FILE \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
| awk '$1=="'$FAMILY'" {print $0}' \
| cat $CORE_PATH/$PROJECT/TEMP/"CONTROL_PED_FILE_FOR_"$FAMILY".ped" /dev/stdin \
>| $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".ped"

#############################################################################
##### 01. Subset BiAllelic SNVs with global MAF from 1000 genomes > 0.1 #####
#############################################################################

### First Reformat the OneKGP.AF tag to OneKGP_AF b/c GATK will not recognize foo.bar correctly and will just look for foo.

sed 's/OneKGP.AF/OneKGP_AF/g' $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.vcf" \
>| $CORE_PATH/$PROJECT/TEMP/VCF_PREP/${FAMILY}.VQSR.PASS.SNV.reformat.vcf

### Extract SNVS with One thousand genome MAF > 0.1 (10 percent)

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-R $REF_GENOME \
-T SelectVariants \
--variant $CORE_PATH/$PROJECT/TEMP/VCF_PREP/${FAMILY}.VQSR.PASS.SNV.reformat.vcf \
--selectexpressions 'OneKGP_AF > 0.1' \
-o $CORE_PATH/$PROJECT/TEMP/VCF_PREP/${FAMILY}.VQSR.PASS.SNV.OneKG_AF.vcf

#########################################
##### 02. Convert vcf to PLINK file #####
#########################################

$VCFTOOLS_DIR/vcftools \
--vcf $CORE_PATH/$PROJECT/TEMP/VCF_PREP/${FAMILY}.VQSR.PASS.SNV.OneKG_AF.vcf \
--plink-tped \
--out $CORE_PATH/$PROJECT/TEMP/PLINK/${FAMILY}.VQSR.PASS.SNV

# vcftools does not know the information here since it is just parsing the vcf file...so Hua is just moving this file out of the way
# The next step after this is really reordering the ped file (here it is called tfam) to match the sample order in the vcf file

mv -v $CORE_PATH/$PROJECT/TEMP/PLINK/${FAMILY}.VQSR.PASS.SNV.tfam \
$CORE_PATH/$PROJECT/TEMP/PLINK/${FAMILY}.VQSR.PASS.SNV.tfam.bak

# Storing the path to concantenated ped file as a variable

CONTROL_PED_FILE_FOR_FAMILY=`echo $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".ped"`

## Replace the tfam file with right pedigree info and in the sample order in the vcf file ##

zgrep -m 1 "^#CHROM" $CORE_PATH/$PROJECT/TEMP/VCF_PREP/"CONTROLS_PLUS_"$FAMILY".VQSR.ANNOTATED.SNV_ONLY.PASS.BIALLELIC.vcf" \
| sed 's/\t/\n/g' \
| awk 'NR>9' \
| awk '{print "awk \x27 $2==\x22"$0"\x22 \x27","'$CONTROL_PED_FILE_FOR_FAMILY'"}' \
| bash \
>| $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.tfam

####################################################
##### 03.A Run Relatedness check using KING1.9 #####
####################################################

##	Pedigree file needs to be modified, a final list of Coriell samples needs to be considered ##

$PLINK2_DIR/plink \
--noweb \
--tfile $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV \
--maf 0.1 \
--make-bed \
--out $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin

$KING_DIR/king -b $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin.bed \
--kinship \
--IBS \
--prefix $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.KinshipIBS

# this is for homogenous populations, but this will never be the case

$KING_DIR/king -b $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin.bed \
 --homo \
 --showIBD \
 --minMAF 0.1 \
 --prefix $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.HomoIBD

awk '$4!=0 {print $0}' $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.HomoIBD.kin \
>| $CORE_PATH/$PROJECT/$FAMILY/RELATEDNESS/$FAMILY.VQSR.PASS.SNV.HomoIBD.kin.final.txt

awk '$4!=0 {print $0}' $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.HomoIBD.kin0 \
>| $CORE_PATH/$PROJECT/$FAMILY/RELATEDNESS/$FAMILY.VQSR.PASS.SNV.HomoIBD.kin0.final.txt

# format and concatenate kinship ibs output

(awk 'BEGIN {print "FID1","ID1","FID2","ID2","Phi","N_IBS0","N_IBS1","N_IBS2","IBS","SE_IBS","N_HetHet","N_Het1","N_Het2","Distance","SE_Dist","Kinship","Error"} \
NR>1 {print $1,$2,$1,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.KinshipIBS.ibs ;
awk 'NR>1 $15>0.06 {print $1,$2,$3,$4,"0",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"1"} \
$15<=0.06 {print $1,$2,$3,$4,"0",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"0"}' $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.KinshipIBS.ibs0) \
| sed 's/ /\t/g' \
>| $CORE_PATH/$PROJECT/$FAMILY/RELATEDNESS/$FAMILY.VQSR.PASS.SNV.KinshipIBS.txt

###########################################
##### Try PLINK for relatedness check #####
###########################################

$PLINK2_DIR/plink --noweb --maf 0.1 --genome --genome-full \
--bfile $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin \
--out $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.genome

############################################################
##### Reformat output to make a better delimited table #####
############################################################

sed -r 's/^ *//g ; s/[[:space:]]+/\t/g' $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.genome.genome \
>| $CORE_PATH/$PROJECT/$FAMILY/RELATEDNESS/$FAMILY.VQSR.PASS.SNV.PLINK2.final.txt

##################################################
##### 03.B Run PCA using KING1.9 on sunrhel4 #####
##################################################

$KING_DIR/king -b $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin.bed \
--mds --ibs \
--prefix $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.MDS_IBS

sed -r 's/^ *//g ; s/[[:space:]]+/\t/g' $CORE_PATH/$PROJECT/TEMP/KING/$FAMILY.VQSR.PASS.SNV.MDS_IBSpc.ped \
>| $CORE_PATH/$PROJECT/$FAMILY/PCA/$FAMILY.VQSR.PASS.SNV.MDS_IBSpc.ped.final.txt

#############################
##### Try PLINK for PCA #####
#############################

$PLINK2_DIR/plink --noweb --bfile $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.bin \
--read-genome $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.genome.genome \
--cluster --mds-plot 4 \
--out $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.mds

############################################################
##### Reformat output to make a better delimited table #####
############################################################

sed -r 's/^ *//g ; s/[[:space:]]+/\t/g' $CORE_PATH/$PROJECT/TEMP/PLINK/$FAMILY.VQSR.PASS.SNV.mds.mds \
>| $CORE_PATH/$PROJECT/$FAMILY/PCA/$FAMILY.VQSR.PASS.SNV.mds.PLINK2.final.txt
