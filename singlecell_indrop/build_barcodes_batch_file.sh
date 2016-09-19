for INPUT_FASTQ_GZ in /data/langenau/singlecell_prkdc/fastq/Zebrafish*.processed.fastq.gz
do
	FNAME=`basename ${INPUT_FASTQ_GZ} | sed 's/\.fastq\.gz//g'`
	BARCODES_ABOVE_THRESHOLD=/data/langenau/singlecell_prkdc/demultiplexed/${FNAME}.barcodes_above_threshold.txt
	
	if [ -f /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/barcodes.batch ]; then
		rm /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/barcodes.batch
	fi
	for valid_barcode in `cut -f1 ${BARCODES_ABOVE_THRESHOLD}`
	do
	        echo ${FNAME}.${valid_barcode}  /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.umi /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${valid_barcode}.fastq.gz  | sed 's/ /\t/g' >> /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/barcodes.batch
	done
done
#rm -rf /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*_fqs
