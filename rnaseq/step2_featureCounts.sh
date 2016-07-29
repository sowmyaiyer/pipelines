ANNOTATION_GTF=$1 #/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf #$1
OUTPUT_DIR=$2 #/data/rivera/lombard/out #$2
STRANDEDNESS=$3 # 0 (unstranded), 1 (stranded) or 2 (reversely stranded)

STAR_out=${OUTPUT_DIR}/STAR_out

mkdir -p  ${OUTPUT_DIR}/featureCounts_out

featureCounts -s ${STRANDEDNESS} -p --primary -Q 10 -T 4 -C -F GTF -a ${ANNOTATION_GTF} -o ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts ${STAR_out}/*Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
sed 1d  ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts | sed 's/\"//g' | sed "s|${STAR_out}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt
sed 's/\"//g' ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.summary | sed "s|${STAR_out}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' >  ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.summary.formatted.txt

sed 1d ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.noheader.txt

head -1 ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.header

cat ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.header ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.noheader.txt > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.txt
echo "done STAR feature counts"

