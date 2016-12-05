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
	#ANNOTATION_GTF=/data/rivera/genomes/mm10/gencode.vM6.basic.annotation.gtf
	ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.mm10.09_30_2016/mm10.09_30_2016.refFlat.gtf
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
elif [[ ${GENOME} == "GRCz10" ]]; then
        GENOME_FASTA=/data/langenau/Danio_rerio.GRCz10.dna.toplevel.fa #/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
        STAR_INDEX_DIR=/data/langenau/STAR_index_GRCz10.85
        ANNOTATION_GTF=/data/langenau/Danio_rerio.GRCz10.85.gtf
        if [[ ! -f  /data/rivera/sowmya/genomes/GRCz10.chrom.sizes ]]; then
                /data/aryee/pub/genomes/fetchChromSizes danRer10 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/GRCz10.chrom.sizes
        fi
fi
echo $GENOME_FASTA $STAR_INDEX_DIR $ANNOTATION_GTF



module use /apps/modulefiles/lab
module load aryee/star-2.4.0h
module load samtools/1.1
module load BEDTools_2.17
module use /apps/modulefiles/test
module load python/2.7.3
module load RSeQC/2.4
module load ucsc
module load cufflinks/2.2.1

reads_after_rRNA_removal=`samtools view -c  ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam`

# Generate depth-normalized bigwig file
echo "start bw file generation"
scalingFactor=`echo ${reads_after_rRNA_removal} | awk -vmappedReads=${reads_after_rRNA_removal} -vSCALING_FACTOR=${SCALING_FACTOR_MULT} 'BEGIN{print 1000000*SCALING_FACTOR/mappedReads}'`
echo $scalingFactor
bedtools genomecov -ibam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam -g  /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes -bga -split -scale ${scalingFactor} >  ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg   /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bw
rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
echo "DONE bw file generation"
