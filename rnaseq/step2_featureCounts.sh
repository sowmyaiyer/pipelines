if [[ $# -lt 2 ]]; then
    echo "Usage: `basename $0` <full path to GTF> <output directory(parent folder to STAR_out from previous step)> <strandedness>(if stranded library, "forward" or "reverse"; if not provided, assumed unstranded)... ">&2
    exit 1
fi

ANNOTATION_GTF=$1 #/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf #$1
OUTPUT_DIR=$2 #/data/rivera/out #$2
STRANDEDNESS=$3 # if stranded library, "forward" or "reverse"; if unstranded, do not provide #3

if [[ ${STRANDEDNESS} == "reverse" ]]; then
	FEATURE_COUNTS_STRANDEDNESS=2
elif [[ ${STRANDEDNESS} == "forward" ]]; then
	FEATURE_COUNTS_STRANDEDNESS=1
else
	FEATURE_COUNTS_STRANDEDNESS=0
fi



STAR_out=${OUTPUT_DIR}/STAR_out

mkdir -p  ${OUTPUT_DIR}/featureCounts_out

featureCounts -s ${FEATURE_COUNTS_STRANDEDNESS} -p --primary -Q 10 -T 4 -C -F GTF -a ${ANNOTATION_GTF} -o ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts ${STAR_out}/*Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
sed 1d  ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts | sed 's/\"//g' | sed "s|${STAR_out}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt
sed 's/\"//g' ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.summary | sed "s|${STAR_out}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' >  ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.summary.formatted.txt

sed 1d ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.noheader.txt

head -1 ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.header

cat ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.header ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.formatted.noheader.txt > ${OUTPUT_DIR}/featureCounts_out/STAR_combined.featureCounts.txt
echo "done STAR feature counts"

