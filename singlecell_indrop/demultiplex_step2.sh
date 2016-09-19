INPUT_FASTQ_GZ=$1
FNAME=`basename ${INPUT_FASTQ_GZ} | sed 's/\.fastq\.gz//g'`
BARCODES_ABOVE_THRESHOLD=/data/langenau/singlecell_prkdc/demultiplexed/${FNAME}.barcodes_above_threshold.txt

mkdir -p /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}
for valid_barcode in `cut -f1 ${BARCODES_ABOVE_THRESHOLD}`
do
        cat /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*_fqs/${valid_barcode}.fq | gzip -c > /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.fastq.gz
        ls  /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.fastq.gz
        cat /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*_fqs/${valid_barcode}.umi > /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.umi
        ls  /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.umi
done
cat /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*_fqs/_barcodes_below_threshold.fq | gzip -c > /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}.barcodes_below_threshold.fastq.gz

#rm -rf /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*_fqs
