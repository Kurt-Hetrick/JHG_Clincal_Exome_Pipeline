# for X it's 1297674, not inside BED
BED=(7821961 5837497 4499324 3091846 3506196 3878112 3658664 2566853 3146790 3101166 4523961 4255702 1432634 2414646 2787139 3189648 4417001 1220227 4581852 1857126 789143 1626513)

# exomeDepth result file 
EXOM=$1
# example command to create bed file
# sed '1d'  RefGene.Coding.Clean.Format.19Feb2018_gene_exomeDepth_noY.bed |sort -k1,1 -k2,2n -k3,3n -u|bedtools merge -i - >RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region.bed
# BEDFILE=/mnt/research/statgen/DDL/CNV_pipeline_DDL/data/RefGene.Coding.Clean.Format.19Feb2018_gene_noY_uniq_region2.bed
BEDFILE=$2
#BEDTOOLS=/cm/shared/apps/clinical/bedtools/2.26.0/bin/bedtools

for i in `seq 22`
do
  # echo "chr "$i 
  # awk -v chr="$i" '{if($1==chr) print $1,$2,$3;OFS="\t" }' ../RefGene.Coding.Clean.Format.19Feb2018_gene_noXY_exomeDepth.bed |sort -k 2,2n |bedtools merge -i - |awk '{SUM += $3-$2} END {print SUM}'
 awk -v chromosome="$i" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' $EXOM |sort -k 2,2n |bedtools merge -i - \
 |bedtools intersect -a - -b $BEDFILE |awk -v len=${BED[$i-1]} -v ch=$i '{SUM += $3-$2} END {printf "chr" ch " " "%.1f\n", 100*SUM/len}'
done

awk -v chromosome="X" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' $EXOM |sort -k 2,2n |bedtools merge -i - |bedtools intersect -a - -b $BEDFILE |awk -v len="2760475" '{SUM += $3-$2} END {printf "chrX " "%.1f\n", 100*SUM/len}'

awk -v chromosome="Y" '{if($1==chromosome) print $1,$2,$3;OFS="\t" }' $EXOM |sort -k 2,2n |bedtools merge -i - |bedtools intersect -a - -b $BEDFILE |awk -v len="158436" '{SUM += $3-$2} END {printf "chrY " "%.1f\n", 100*SUM/len}'