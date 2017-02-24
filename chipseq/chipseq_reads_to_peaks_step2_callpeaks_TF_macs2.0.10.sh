if [[ $# -lt 4 ]]; then
    echo "Usage: `basename $0` 
	-s=sampleName 
	-t=<treatment bam file> 
	-c=<control bam file> 
	-q=<q-value threshold, 0.01 default if not provided> 
	-o=<outputDir> 
	-d=(merge peaks within d; optional) 
	--broad(optional) ... ">&2
    exit 1
fi

QTHRESH=0.01
BROAD_PEAKS="NO"
MERGE_PEAKS="NO"

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
    -d=*|--mergeDistance=*)
    MERGE_PEAKS="YES"
    MERGE_DISTANCE="${i#*=}"
    shift # past argument with no value
    ;;
   -h)
    echo "Usage: `basename $0`
        -s=sampleName
        -t=<treatment bam file>
        -c=<control bam file>
        -q=<q-value threshold, 0.01 default if not provided>
        -o=<outputDir>
        -d=(merge peaks within d; optional)
        --broad(optional) ... ">&2
    exit 1
    shift
    ;;
    *)
    ;;
esac
done


module use /apps/modulefiles/test
module load python/2.7.3
module load  numpy-1.7
module load macs2/2.0.10

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

if [[ ${BROAD_PEAKS} == "YES" ]]; then
	macs_command="macs2 callpeak --nomodel -g hs -q ${QTHRESH} -t ${TREATMENT_FILE} -c ${CONTROL_FILE} -n ${SAMPLE_NAME} --broad"
	PEAK_FILE=${OUTPUT_DIR}/${SAMPLE_NAME}_broad_peaks.bed
else
	macs_command="macs2 callpeak -g hs -q ${QTHRESH} -t ${TREATMENT_FILE} -c ${CONTROL_FILE} -n ${SAMPLE_NAME}"
	PEAK_FILE=${OUTPUT_DIR}/${SAMPLE_NAME}_peaks.bed
fi

echo $macs_command

${macs_command} 
PEAK_FILE_PROPER=${PEAK_FILE}.clean.bed
BLACKLIST=/data/rivera/genomes/wgEncodeDacMapabilityConsensusExcludable.bed

# Get valid chromosomes and subtract blacklisted regions
awk 'NR==FNR{a[$0];next}($1 in a)' /PHShome/si992/commonscripts/chipseq/human_chrs.nochr.txt ${PEAK_FILE} | bedtools subtract -a stdin -b ${BLACKLIST} -A > ${PEAK_FILE_PROPER}
echo "DONE cleaning up peaks"

if [[ ${MERGE_PEAKS} == "YES" ]]; then
	sort -k1,1 -k2,2n ${PEAK_FILE_PROPER} | bedtools merge -i stdin -d ${MERGE_DISTANCE} -nms > ${PEAK_FILE}.clean.merged.${MERGE_DISTANCE}.bed
	echo "DONE merging peaks"
fi
