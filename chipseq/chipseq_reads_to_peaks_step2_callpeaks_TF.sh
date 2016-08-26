if [[ $# -lt 4 ]]; then
    echo "Usage: `basename $0` -s=sampleName -t=<treatment bam file> -c=<control bam file> -q=<q-value threshold, 0.01 default> -o=<outputDir> --broad(optional) ... ">&2
    exit 1
fi

QTHRESH=0.01
BROAD_PEAKS="NO"

for i in "$@"
do
case $i in
    -s=*|--sampleName=*)
    SAMPLE_NAME="${i#*=}"
    shift
    ;;
    -t=*|--treatment_bam=*)
    TREATMENT_FILE="${i#*=}"
    shift
    ;;
    -c=*|--control_bam=*)
    CONTROL_FILE="${i#*=}"
    shift
    ;;
    -q=*|--qthresh=*)
    QTHRESH="${i#*=}"
    shift
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"
    shift
    ;;
    --broad)
    BROAD_PEAKS="YES"
    shift # past argument with no value
    ;;
   -h)
    echo "Usage: `basename $0` -s=sampleName -t=<treatment bam file> -c=<control bam file> -q=<q-value threshold> -o=<outputDir> --broad(optional) ... ">&2
    exit 1
    shift
    ;;
    *)
    ;;
esac
done


module load python/2.7.3
module load macs2/2.1

if [[ $BROAD_PEAKS == "YES" ]]; then
	macs_command="macs2 callpeak --nomodel -g hs -q ${QTHRESH} -t ${TREATMENT_FILE} -c ${CONTROL_FILE} -n ${SAMPLE_NAME} --outdir ${OUTPUT_DIR} --broad"
	PEAK_FILE=${OUTPUT_DIR}/${SAMPLE_NAME}_peaks.narrowPeak
else
	macs_command="macs2 callpeak -g hs -q ${QTHRESH} -t ${TREATMENT_FILE} -c ${CONTROL_FILE} -n ${SAMPLE_NAME} --outdir ${OUTPUT_DIR}"
	PEAK_FILE=${OUTPUT_DIR}/${SAMPLE_NAME}_peaks.broadPeak
fi
echo $macs_command

${macs_command} 
PEAK_FILE_PROPER=${PEAK_FILE}.clean.bed
BLACKLIST=/data/rivera/genomes/wgEncodeDacMapabilityConsensusExcludable.bed

# Get valid chromosomes and subtract blacklisted regions
awk 'NR==FNR{a[$0];next}($1 in a)' $HOME/commonscripts/chipseq/human_chrs.nochr.txt ${OUTPUT_DIR}/${SAMPLE_NAME}_peaks.narrowPeak | bedtools subtract -a stdin -b ${BLACKLIST} -A > ${OUTPUT_DIR}/${SAMPLE_NAME}_peaks.narrowPeak.clean.bed

echo "DONE cleaning up peaks"
