#GENOME=$1 #Zv9
#FILE_SAMPLE_LIST=$2 #test_Zv9_config.txt
#OUTPUT_DIR=$3 # /data/langenau/out
#COMBINE_FASTQ=$4 # true or false


if [[ $# -lt 3 ]]; then
    echo "Usage: `basename $0` -g=genome -s=sampleFile -o=outputDir -r(if stranded library, FR or RF; if unstranded, do not provide) --combineFastq(if lanes are to be combined)... ">&2
    exit 1
fi

COMBINE_FASTQ="NO"
STRANDED="NO"

for i in "$@"
do
case $i in
    -g=*|--genome=*)
    GENOME="${i#*=}" #Zv9
    shift # past argument=value
    ;;
    -s=*|--sampleFile=*)
    FILE_SAMPLE_LIST="${i#*=}" #test_Zv9_config.txt
    shift # past argument=value
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"  # /data/langenau/out
    shift # past argument=value
    ;;
    -r=*|--strandedness=*)
    STRANDED="YES"
    STRANDEDNESS="${i#*=}"  # for hisat2  --rna-strandness FR or RF. Default unstranded
    shift # past argument=value
    ;;
    --combineFastq)
    COMBINE_FASTQ="YES"
    shift # past argument with no value
    ;;
    *)
            # unknown option
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
elif [[ ${GENOME} == "Zv9" ]]; then
	GENOME_FASTA=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
	STAR_INDEX_DIR=/data/langenau/STAR_index_Zv9
        ANNOTATION_GTF=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf 
	if [[ ! -f  /data/rivera/sowmya/genomes/Zv9.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes danRer7 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/Zv9.chrom.sizes
	fi
fi
echo $GENOME_FASTA $STAR_INDEX_DIR $ANNOTATION_GTF

while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`

	echo """
	module load aryee/star-2.4.0h
	module load samtools/1.1

	module load python/2.7.3
        module load RSeQC/2.4
        module load ucsc
	module load cufflinks/2.2.1

	# Remove rRNA reads

        split_bam.py -i ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam -r ${HOME}/commonscripts/rnaseq/${GENOME}_rRNA.bed --out-prefix=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed > ${OUTPUT_DIR}/STAR_out/rRNA_${SAMPLE_NAME}.txt
        mv ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.ex.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
        samtools index ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
        rm ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.in.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.junk.bam

	# Get some numbers on mapping

	totalRawReads=\`grep \"Number of input reads\" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2\`
	uniquelyMappedReads=\`grep \"Uniquely mapped reads %\" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2\`
	multiMappingReads=\`grep \"% of reads mapped to multiple loci\" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Log.final.out | cut -f2\`
        dedupedReads=\`samtools view -c  ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam\`
        fraction_dups=\`grep \"Unknown \" ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.dups.txt | awk '{ print \$(NF-1)}'\`
        reads_after_rRNA_removal=\`samtools view -c  ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam \`
        echo ${SAMPLE_NAME} \${totalRawReads} \${uniquelyMappedReads} \${multiMappingReads} \${dedupedReads} \${fraction_dups} \${reads_after_rRNA_removal} > ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.mapping_metrics.txt


	# Generate depth-normalized bigwig file

        scalingFactor=\`echo \${mappedReads} | awk -vmappedReads=\${mappedReads} 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam -g  /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes -bga -split -scale ${scalingFactor} >  ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
        bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg   /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bw
	rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg

	# Get FPKMs for isoforms
	cufflinks -o ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME} -G ${ANNOTATION_GTF} ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
	
	"""  > ../bsubFiles/star_alignment_${GENOME}_${SAMPLE_NAME}.bsub
done < ${FILE_SAMPLE_LIST}


