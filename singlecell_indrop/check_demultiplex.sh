INPUT_FASTQ_GZ=$1
FNAME=`basename ${INPUT_FASTQ_GZ} | sed 's/\.fastq\.gz//g'`
CB_HISTOGRAM_FILE=`basename ${INPUT_FASTQ_GZ} | sed 's/\.processed\.fastq\.gz/.cb_histogram.txt/g' `
THRESHOLD=10000

for lineall in `awk '{ if ($2 > 10000) print $1"_"$2}' /data/langenau/singlecell_prkdc/fastq/$CB_HISTOGRAM_FILE`
do
	barcode=`echo "$lineall" | awk -F"_" '{ print $1}'`
	linecount=`echo "$lineall" | awk -F"_" '{ print 4*$2}'`
	fastq_linecount=`zcat /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/${barcode}.fastq.gz | wc -l`
	if [[ $linecount -ne $fastq_linecount ]]; then
		echo uh-oh! $barcode numbers not matching  $linecount $fastq_linecount
	fi
done 
