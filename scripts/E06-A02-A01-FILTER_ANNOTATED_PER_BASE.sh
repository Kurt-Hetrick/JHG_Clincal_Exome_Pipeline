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

	CORE_PATH=$1
	PROJECT=$2
	FAMILY=$3
	SM_TAG=$4
	CODING_BED=$5
		CODING_BED_NAME=$(basename ${CODING_BED} .bed)
		CODING_MD5=$(md5sum ${CODING_BED} | cut -c 1-7)
	PADDING_LENGTH=$6

# The input bed file could be a variable name based on the padding length
# Remove the bases in the annotated per base report to 

START_PER_BASE_FILTER=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# bases below 30x
# uncommitted to this being a variable, but it can change

	(head -n 1 ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt ; \
	awk 'NR>1' ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.txt \
		| awk '$7<30' ) \
	>| ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/DEPTH_OF_COVERAGE/CODING_PADDED/${SM_TAG}_${CODING_BED_NAME}-${CODING_MD5}-${PADDING_LENGTH}-BP-PAD.PER.BASE.REPORT.lt30.txt

END_PER_BASE_FILTER=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${SM_TAG}_${PROJECT},G01,REFSEQ_PER_BASE_FILTER,${HOSTNAME},${START_PER_BASE_FILTER},${END_PER_BASE_FILTER} \
	>> ${CORE_PATH}/${PROJECT}/${FAMILY}/${SM_TAG}/REPORTS/${PROJECT}.WALL.CLOCK.TIMES.csv