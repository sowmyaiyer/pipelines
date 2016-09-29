if [[ $# -lt 4 ]]; then
    echo "Usage: `basename $0` 
	-g=genome 
	-s=sampleName 
	-q=FastqPrefix 
	-o=outputDir 
	-r=library strandedness(optional; forward or reverse; defaults to unstranded)
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
    -r=*|--strandedness=*)
    STRANDEDNESS="${i#*=}"  
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
        -r=library strandedness(optional; forward or reverse; defaults to unstranded)
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
	STAR_INDEX_DIR=/data/rivera/genomes/mm10/star_index
	ANNOTATION_GTF=/data/rivera/genomes/mm10/gencode.vM6.basic.annotation.gtf
	if [[ ! -f  /data/rivera/sowmya/genomes/mm10.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes mm10  > /data/rivera/sowmya/genomes/mm10.chrom.sizes
	fi
elif [[ ${GENOME} == "hg19" ]]; then
	GENOME_FASTA=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
	STAR_INDEX_DIR=/data/rivera/sowmya/genomes/STAR_indices/STAR_index_hg19
        ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf
	if [[ ! -f  /data/rivera/sowmya/genomes/hg19.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes hg19 > /data/rivera/sowmya/genomes/hg19.chrom.sizes
	fi
fi
echo $GENOME_FASTA $STAR_INDEX_DIR $ANNOTATION_GTF

if [[ ${STRANDEDNESS} == "reverse" ]]; then
        CUFFLINKS_STRANDEDNESS="fr-firststrand"
elif [[ ${STRANDEDNESS} == "forward" ]]; then
        CUFFLINKS_STRANDEDNESS="fr-secondstrand"
else
        CUFFLINKS_STRANDEDNESS="fr-unstranded"
fi

mkdir -p ${OUTPUT_DIR}/STAR_out
mkdir -p ${OUTPUT_DIR}/cufflinks_out
mkdir -p ${OUTPUT_DIR}/bigwigs


R1_FASTQ=`ls ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz`
R2_FASTQ=`ls ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz`

echo sample name is ${SAMPLE_NAME}
echo fastq_1 is ${R1_FASTQ}
echo fastq_2 is ${R2_FASTQ}

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

echo "Start STAR alignment"
readLength=`zcat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz | head -2 | tail -1 | awk '{ print length($0)}'` # only checks read length of first read. Maybe a problem if reads have been trimmed or have variable lengths for other reasons
echo read length is $readLength
sjdbOverhang=$((readLength -1))
star_fastq_R1=`echo $R1_FASTQ | sed 's/ /,/g'`
star_fastq_R2=`echo $R2_FASTQ | sed 's/ /,/g'`

module use /apps/modulefiles/lab
module load aryee/star-2.4.0h
module load samtools/1.1
module load BEDTools_2.17
module use /apps/modulefiles/test
module load python/2.7.3
module load RSeQC/2.4
module load ucsc
module load cufflinks/2.2.1

star_fastq_R1=`echo $R1_FASTQ | sed 's/ /,/g'`
star_fastq_R2=`echo $R2_FASTQ | sed 's/ /,/g'`	

if [[ $CUFFLINKS_STRANDEDNESS == "fr-unstranded" ]]; then
	added_star_option=" --outSAMstrandField intronMotif"
else
	added_star_option=""
fi

STAR ${added_star_option} --genomeDir $STAR_INDEX_DIR --genomeFastaFiles $GENOME_FASTA --sjdbGTFfile $ANNOTATION_GTF --sjdbOverhang ${sjdbOverhang} --runThreadN 1 --outSAMtype BAM Unsorted --outFileNamePrefix  ${OUTPUT_DIR}/STAR_out/$SAMPLE_NAME --readFilesCommand zcat --readFilesIn $star_fastq_R1 $star_fastq_R2

echo "Done STAR alignment"
	
# Sorting in separate step because STAR sometimes crashes while sorting by coordinate
echo "Start Sorting BAM file"

samtools sort -T ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out -o ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.out.bam

echo "End Sorting BAM file"

# Remove duplicates
echo "Start deduping"

java  -Xmx40g -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam O=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam M=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR

echo "DONE deduping"


# Remove rRNA reads
echo "start removing rRNA reads"
split_bam.py -i ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam -r /PHShome/si992/commonscripts/rnaseq/${GENOME}_rRNA.bed --out-prefix=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed > ${OUTPUT_DIR}/STAR_out/rRNA_${SAMPLE_NAME}.txt
mv ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.ex.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
samtools index ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
rm ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.in.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.junk.bam
echo "DONE removing rRNA reads"

# Get some numbers on mapping
echo "start getting mapping metrics"
totalRawReads=`grep "Number of input reads" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2`
uniquelyMappedReads=`grep "Uniquely mapped reads %" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2`
multiMappingReads=`grep "% of reads mapped to multiple loci" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2`
dedupedReads=`samtools view -c  ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam`
fraction_dups=`grep "Unknown " ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.dups.txt | awk '{ print $(NF-1)}'`
reads_after_rRNA_removal=`samtools view -c  ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam `
echo ${SAMPLE_NAME} ${totalRawReads} ${uniquelyMappedReads} ${multiMappingReads} ${dedupedReads} ${fraction_dups} ${reads_after_rRNA_removal} > ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.mapping_metrics.txt
echo "DONE mapping metrics"

# Generate depth-normalized bigwig file
echo "start bw file generation"
scalingFactor=`echo ${reads_after_rRNA_removal} | awk -vmappedReads=${reads_after_rRNA_removal} -vSCALING_FACTOR=${SCALING_FACTOR_MULT} 'BEGIN{print 1000000*SCALING_FACTOR/mappedReads}'`
echo $scalingFactor
bedtools genomecov -ibam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam -g  /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes -bga -split -scale ${scalingFactor} >  ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg   /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bw
rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
echo "DONE bw file generation"

# Get FPKMs for isoforms
echo "start cufflinks"
cufflinks --quiet --library-type ${CUFFLINKS_STRANDEDNESS} -o ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME} -G ${ANNOTATION_GTF} ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
echo "DONE cufflinks"
# DONE
echo "ALL DONE!" See ${OUTPUT_DIR} for output and logs
