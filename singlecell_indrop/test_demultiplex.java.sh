INPUT_FASTQ_GZ=$1
CB_FILE=$2
OUTPUT_DIR=$3

if [ -z "$INPUT_FASTQ_GZ" ]; then
        echo "No fastq.gz file supplied. Exiting..."
        exit 1
fi
if [ ! -f "$INPUT_FASTQ_GZ" ]; then
        echo "Fastq file $INPUT_FASTQ_GZ does not exist. Exiting..."
        exit 1
fi

if [ -z "$CB_FILE" ]; then
        echo "No barcode histogram file name supplied, Exiting..."
        exit 1
fi
if [ ! -f "$CB_FILE" ]; then
        echo "barcode histogram file $CB_FILE does not exist, Exiting.. "
        exit 1
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo "output directory not provided, Exiting.. "
    exit 1
fi


THRESHOLD=10000
FNAME=`basename ${INPUT_FASTQ_GZ} | sed 's/\.fastq\.gz//g'`

mkdir -p /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/

# Following line needs to be executed only once
#zcat ${INPUT_FASTQ_GZ} | split - --suffix-length=5 --lines=1000000 --numeric-suffixes /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/
echo "done splitting"

awk -vthresh=${THRESHOLD} '{ if ($2 > thresh) print $1 }' ${CB_FILE} > ${OUTPUT_DIR}/${FNAME}.barcodes_above_threshold.txt
for split_file in /data/langenau/singlecell_prkdc/processing/demultiplex_splits4/${FNAME}.split/*[0-9]
do
	echo ${split_file}
	mkdir -p ${split_file}_fqs
	echo """
	java -cp /PHShome/si992/commonscripts/singlecell_indrop/samtools-htsjdk-c3d5a88-unspecified-SNAPSHOT.jar:/PHShome/si992/commonscripts/singlecell_indrop demultiplex.TryFastqReader ${split_file} ${OUTPUT_DIR}/${FNAME}.barcodes_above_threshold.txt ${split_file}_fqs/
	""" > ../bsubFiles/demultiplex.${FNAME}.`basename ${split_file}`.bsub
done 
