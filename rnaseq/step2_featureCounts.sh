if [[ $# -lt 2 ]]; then
	echo "Usage: `basename $0` 
	-g=<genome>
	-b=<txt file with list of bams> NOT IMPLEMENTED
	-p=<path to bam files for wildcard search>
	-o=<featureCounts output folder>
	-r=<strandedness>(if stranded library, "forward" or "reverse"; if not provided, assumed unstranded)... ">&2    
	exit 1
fi


for i in "$@"
do
case $i in
    -g=*|--genome=*)
    GENOME="${i#*=}"
    shift
    ;;
    -b=*|--bamFileList=*)
    BAMFILE_LIST="${i#*=}"
    shift
    ;;
    -p=*|--bamFilePath=*)
    BAMFILE_PATH="${i#*=}"
    shift
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"
    shift
    ;;
    -r=*|--strandedness=*)
    STRANDEDNESS="${i#*=}"
    shift
    ;;
    -h)
        "Usage: `basename $0`
        -g=<genome>
        -b=<txt file with list of bams> NOT IMPLEMENTED
        -p=<path to bam files for wildcard search>
        -o=<featureCounts output folder>
        -r=<strandedness>(if stranded library, "forward" or "reverse"; if not provided, assumed unstranded)... ">&2
	exit 1
    shift
    ;;
    *)
    ;;
esac
done


if [[ ${GENOME} == "mm10" ]]; then
        ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.mm10.09_30_2016/mm10.09_30_2016.refFlat.gtf
elif [[ ${GENOME} == "mm9" ]]; then
        ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.mm9.12_05_2016/mm9.12_05_2016.refFlat.gtf
elif [[ ${GENOME} == "hg19" ]]; then
        ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf
elif [[ ${GENOME} == "GRCz10" ]]; then
        ANNOTATION_GTF=/data/langenau/Danio_rerio.GRCz10.85.gtf
fi


if [[ ${STRANDEDNESS} == "reverse" ]]; then
	FEATURE_COUNTS_STRANDEDNESS=2
elif [[ ${STRANDEDNESS} == "forward" ]]; then
	FEATURE_COUNTS_STRANDEDNESS=1
else
	FEATURE_COUNTS_STRANDEDNESS=0
fi

mkdir -p ${OUTPUT_DIR}

featureCounts -s ${FEATURE_COUNTS_STRANDEDNESS} -p --primary -Q 10 -T 4 -C -F GTF -a ${ANNOTATION_GTF} -o ${OUTPUT_DIR}/combined.counts ${BAMFILE_PATH}/*Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
sed 1d  ${OUTPUT_DIR}/combined.counts | sed 's/\"//g' | sed "s|${BAMFILE_PATH}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' > ${OUTPUT_DIR}/combined.counts.formatted.txt
sed 's/\"//g' ${OUTPUT_DIR}/combined.counts.summary | sed "s|${BAMFILE_PATH}\/*||g"| sed 's/Aligned\.sortedByCoord\.out\.deduped\.rRNA_removed\.bam//g' >  ${OUTPUT_DIR}/combined.counts.summary.formatted.txt

sed 1d ${OUTPUT_DIR}/combined.counts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/combined.counts.formatted.noheader.txt

head -1 ${OUTPUT_DIR}/combined.counts.formatted.txt | cut -f1,7- > ${OUTPUT_DIR}/combined.counts.header

cat ${OUTPUT_DIR}/combined.counts.header ${OUTPUT_DIR}/combined.counts.formatted.noheader.txt > ${OUTPUT_DIR}/combined.counts.txt

rm ${OUTPUT_DIR}/combined.counts.formatted.txt ${OUTPUT_DIR}/combined.counts.formatted.noheader.txt ${OUTPUT_DIR}/combined.counts.header ${OUTPUT_DIR}/combined.counts.summary
echo "done feature counts"

