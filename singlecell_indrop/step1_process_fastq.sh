if [[ $# -lt 2 ]]; then
    echo "Usage: `basename $0` <sample file> <output dir> ...">&2
    exit 1
fi

SAMPLE_FILE=$1
OUTPUT_DIR=$2
TRANSFORM_FILE=$HOME/commonscripts/singlecell_indrop/harvard-indrop-v2.json # Make this $3 later to accomodate for changing protocols
while read line_all
do
	SAMPLE_NAME=`echo $line_all | awk '{ print $1}' `
	FASTQ_READ1=`echo $line_all | awk '{ print $2}' `
	FASTQ_READ2=`echo $line_all | awk '{ print $3}' `
	echo """
	umis fastqtransform --dual_index --separate_cb ${TRANSFORM_FILE} ${FASTQ_READ1}  ${FASTQ_READ2} > ${OUTPUT_DIR}/${SAMPLE_NAME}.parsed.fastq
	umis cb_filter --nedit 1 --bc1 /data/langenau/singlecell_prkdc/gel_barcode1_list_revComp.txt --bc2 /data/langenau/singlecell_prkdc/gel_barcode2_list_revComp.txt ${OUTPUT_DIR}/${SAMPLE_NAME}.parsed.fastq > ${OUTPUT_DIR}/${SAMPLE_NAME}.processed.fastq
	umis cb_histogram ${OUTPUT_DIR}/${SAMPLE_NAME}.processed.fastq > ${OUTPUT_DIR}/${SAMPLE_NAME}.cb_histogram.txt	
	umis cb_histogram ${OUTPUT_DIR}/${SAMPLE_NAME}.parsed.fastq > ${OUTPUT_DIR}/${SAMPLE_NAME}.cb_histogram_unfiltered.txt	
	umis umi_histogram ${OUTPUT_DIR}/${SAMPLE_NAME}.processed.fastq > ${OUTPUT_DIR}/${SAMPLE_NAME}.umi_histogram.txt

	gzip ${OUTPUT_DIR}/${SAMPLE_NAME}.processed.fastq
	""" > ../bsubFiles/process_indrop_fastq.${SAMPLE_NAME}.bsub
done < ${SAMPLE_FILE}
