if [[ $# -lt 4 ]]; then
    echo "Usage: `basename $0`
        -g=genome
        -s=sampleName
        -q=FastqPrefix
        -o=outputDir
        -n=scale for bw(optional; multiple of 1 million; defaults to 1)
        --run_fastqc(optional) ... ">&2
    exit 1
fi

FASTQC="NO"
SCALING_FACTOR_MULT=1

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
    -n=*|--scaling_factor=*)
    SCALING_FACTOR_MULT="${i#*=}"
    shift
    ;;
    --run_fastqc)
    FASTQC="YES"
    shift # past argument with no value
    ;;
   -h)
        echo "Usage: `basename $0`
        -g=genome
        -s=sampleName
        -q=FastqPrefix
        -o=outputDir
        -n=scale for bw(optional; multiple of 1 million; defaults to 1)
        --run_fastqc(optional) ... ">&2
    exit 1
    shift
    ;;
    *)
    ;;
esac
done

if [[ ${GENOME} == "mm10" ]]; then
	GENOME_FASTA=/data/rivera/genomes/mm10/mm10.fa
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/mm10
	if [[ ! -f  /data/rivera/sowmya/genomes/mm10.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes mm10  > /data/rivera/sowmya/genomes/mm10.chrom.sizes
	fi
elif [[ ${GENOME} == "hg19" ]]; then
	GENOME_FASTA=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
	BWA_INDEX=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
	if [[ ! -f  /data/rivera/sowmya/genomes/hg19.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes hg19 > /data/rivera/sowmya/genomes/hg19.chrom.sizes
	fi
elif [[ ${GENOME} == "GRCh37" ]]; then
        GENOME_FASTA=/pub/genome_references/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
        BWA_INDEX=/pub/genome_references/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
        if [[ ! -f  /data/rivera/sowmya/genomes/GRCh37.chrom.sizes ]]; then
                cut -f1,2 /pub/genome_references/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa.fai > /data/rivera/sowmya/genomes/GRCh37.chrom.sizes
        fi
elif [[ ${GENOME} == "Zv9" ]]; then
	GENOME_FASTA=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/Zv9
	if [[ ! -f  /data/rivera/sowmya/genomes/Zv9.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes danRer7 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/Zv9.chrom.sizes
	fi
fi


mkdir -p ${OUTPUT_DIR}/fastqc_out
mkdir -p ${OUTPUT_DIR}/bwa_out
mkdir -p ${OUTPUT_DIR}/bigwigs
mkdir -p ${OUTPUT_DIR}/combined_fastqs

if [[ $FASTQC == "YES" ]]; then
	module load fastqc
	mkdir -p ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1
	mkdir -p ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2
	zcat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz | fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1 stdin
	zcat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz | fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2 stdin
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1/stdin_fastqc.zip ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R1_fastqc.zip
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2/stdin_fastqc.zip ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R2_fastqc.zip
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R1/stdin_fastqc.html ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R1_fastqc.html
	mv ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R2/stdin_fastqc.html ${OUTPUT_DIR}/fastqc_out/${SAMPLE_NAME}_R2_fastqc.html
	rm -rf ${OUTPUT_DIR}/fastqc_out_${SAMPLE_NAME}_R[1,2]
		
fi

module load bwa
module load ucsc

cat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
cat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz

R1_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
R2_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz

MIN_QUAL=5

echo "bwa aln -t 4 ${BWA_INDEX} ${R1_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai"
bwa aln -t 4 ${BWA_INDEX} ${R1_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai
bwa aln -t 4 ${BWA_INDEX} ${R2_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai
bwa sampe ${BWA_INDEX} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${R1_FASTQ} ${R2_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam
samtools view -bS ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam
samtools view -b -F2304 -q ${MIN_QUAL} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
samtools sort -T ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted -o ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
java  -Xmx40g -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam O=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam M=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR

samtools index ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam
rm ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam 

samtools flagstat ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam >  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat
totalRawReads=`awk '{ if (NR == 1) print $1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat`
mappedReads=`awk '{ if (NR == 5) print $1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat`
dedupedReads=`samtools view -c  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam`
fraction_dups=`grep "Unknown " ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt | awk '{ print $(NF-1)}'`

rm ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam

echo ${SAMPLE_NAME} ${totalRawReads} ${mappedReads} ${dedupedReads} ${fraction_dups} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.mapping_metrics.txt

echo "start bw file generation"
scalingFactor=`echo ${dedupedReads} | awk -vmappedReads=${dedupedReads} -vSCALING_FACTOR=${SCALING_FACTOR_MULT} 'BEGIN{print 1000000*SCALING_FACTOR/mappedReads}'`
echo $scalingFactor
bedtools bamtobed -i ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam | bedtools genomecov -i stdin -bga -scale ${scalingFactor} -g /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes >  ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg   /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bw
rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
echo "DONE bw file generation"
