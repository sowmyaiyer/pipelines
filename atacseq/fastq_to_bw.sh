if [[ $# -lt 4 ]]; then
    echo "Usage: `basename $0` 
	-g=genome 
	-s=sampleName 
	-q=FastqPrefix 
	-o=outputDir 
	--combine_fastq(optional)
	--run_fastqc(optional) ... ">&2
    exit 1
fi

FASTQC="NO"
COMBINE_FASTQ="NO"

for i in "$@"
do
case $i in
    -g=*|--genome=*)
    GENOME="${i#*=}" 
    shift 
    ;;
    -s=*|--sampleName=*)
    SAMPLE_NAME="${i#*=}"
    shift 
    ;;
    -q=*|--fastq_prefix=*)
    FASTQ_PREFIX="${i#*=}"
    shift 
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"  
    shift 
    ;;
    --run_fastqc)
    FASTQC="YES"
    shift # past argument with no value
    ;;
    --combine_fastq)
    COMBINE_FASTQ="YES"
    shift # past argument with no value
    ;;
   -h)
	echo "Usage: `basename $0`
        -g=genome
        -s=sampleName
        -q=FastqPrefix
        -o=outputDir
        --run_fastqc(optional) ... ">&2
    exit 1
    shift 
    ;;
    *)
    ;;
esac
done


if [[ ${GENOME} == "mm10" ]]; then
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/mm10
	if [[ ! -f  /data/rivera/sowmya/genomes/mm10.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes mm10  > /data/rivera/sowmya/genomes/mm10.chrom.sizes
	fi
elif [[ ${GENOME} == "hg19" ]]; then
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/hg19
	if [[ ! -f  /data/rivera/sowmya/genomes/hg19.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes hg19 > /data/rivera/sowmya/genomes/hg19.chrom.sizes
	fi
fi
echo $GENOME_FASTA 

mkdir -p ${OUTPUT_DIR}/bwa_out
mkdir -p ${OUTPUT_DIR}/bigwigs

MIN_QUAL=5


echo sample name is ${SAMPLE_NAME}

if [[ $FASTQC == "YES" ]]; then
	echo "start fastqc"
	mkdir -p ${OUTPUT_DIR}/fastqc_out
	module load fastqc/0.11.2
	mkdir -p ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1
	mkdir -p ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2
	zcat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz | fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1 stdin
	zcat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz | fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2 stdin
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1/stdin_fastqc.zip ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R1_fastqc.zip
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2/stdin_fastqc.zip ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R2_fastqc.zip
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1/stdin_fastqc.html ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R1_fastqc.html
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2/stdin_fastqc.html ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R2_fastqc.html
	rm -rf ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R[1,2]
	echo "end fastqc"
fi

echo "Start bwa alignment"

module use /apps/modulefiles/lab
module load aryee/bwa-0.7.12
module load samtools/1.1
module load BEDTools_2.17
module use /apps/modulefiles/test
module load python/2.7.3
module load RSeQC/2.4
module load ucsc

if [[ $COMBINE_FASTQ == "YES" ]]; then
	mkdir -p ${OUTPUT_DIR}/combined_fastq
	cat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz >  ${OUTPUT_DIR}/combined_fastq/${SAMPLE_NAME}.R1.fastq.gz
	cat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz >  ${OUTPUT_DIR}/combined_fastq/${SAMPLE_NAME}.R2.fastq.gz

	R1_FASTQ=${OUTPUT_DIR}/combined_fastq/${SAMPLE_NAME}.R1.fastq.gz
	R2_FASTQ=${OUTPUT_DIR}/combined_fastq/${SAMPLE_NAME}.R2.fastq.gz
else
	R1_FASTQ=${FASTQ_PREFIX}.R1.fastq.gz
	R2_FASTQ=${FASTQ_PREFIX}.R2.fastq.gz
fi

ls ${R1_FASTQ} ${R2_FASTQ}
bwa aln -t 4 ${BWA_INDEX} ${R1_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai
bwa aln -t 4 ${BWA_INDEX} ${R2_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai
bwa sampe ${BWA_INDEX} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${R1_FASTQ} ${R2_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam
samtools view -bS ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam
echo "Done bwa alignment"
samtools view -b -f2 -F2304 -q ${MIN_QUAL} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
echo "start sorting"
samtools sort -T ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted -o ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
echo "done sorting"


echo "Start deduping; index deduped bam"
java  -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam O=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.nodups.bam M=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR
samtools index ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.nodups.bam
echo "done deduping"


echo "Start plotting insert size histogram"
java -jar /apps/source/picardtools/picard-tools-1.58/CollectInsertSizeMetrics.jar I=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.nodups.bam  O=${OUTPUT_DIR}/bwa_out/InsertSizes.${SAMPLE_NAME}.txt H=${OUTPUT_DIR}/bwa_out/insert_size_histogram.${SAMPLE_NAME}.pdf
echo "Done plotting insert size histogram"


echo "start mapping metrics"
samtools flagstat ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam >  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat
totalRawReads=`awk '{ if (NR == 1) print \$1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat`
mappedReads=`awk '{ if (NR == 5) print \$1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat`
dedupedReads=`samtools view -c  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.nodups.bam`
fraction_dups=`grep "Unknown " ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt | awk '{ print $(NF-1)}'`

echo ${SAMPLE_NAME} ${totalRawReads} ${mappedReads} ${dedupedReads} ${fraction_dups} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.mapping_metrics.txt
echo "done mapping metrics"
#rm ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam

echo "start shift reads by +4 and -5 for positive and negative strand reads"
bedtools bamtobed -i ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.nodups.bam | awk 'BEGIN{OFS="\t"} { 
									if ($6 == "+") 
										pos = $2 + 4
									else if ($6 == "-")
										pos = $3 - 5
									print $1,pos,pos+1
									}' > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.shifted.bed
echo "done shift reads"

	
# Generate depth-normalized bigwig file
echo "start bw file generation"
normFactor=`echo "10000000/${dedupedReads}" | bc -l`
echo $normFactor
bedtools coverage -a ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.shifted.bed -b /data/rivera/genomes/hg19.windows.150.20.chr.bed -counts > ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.raw.bg
awk -vnormFactor=${normFactor} 'BEGIN{OFS="\t"}{if ($4 > 0 && $3-$2 > 140) print $1,$2+65,$3-65,normFactor*$4}' ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.raw.bg > ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.unsorted.bg
sort -k1,1 -k2,2n ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.unsorted.bg > ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.bg
bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.bg  /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.bw
#rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.raw.bg ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.unsorted.bg ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.bg
echo "DONE bw file generation"

# DONE
echo "ALL DONE!" See ${OUTPUT_DIR} for output and logs
