#GENOME=$1 #Zv9
#FILE_SAMPLE_LIST=$2 #test_Zv9_config.txt
#OUTPUT_DIR=$3 # /data/langenau/out
#COMBINE_FASTQ=$4 # true or false


if [[ $# -lt 3 ]]; then
    echo "Usage: `basename $0` -g=genome -s=sampleFile -o=outputDir ... ">&2
    exit 1
fi

COMBINE_FASTQ="NO"
STRANDED="NO"

for i in "$@"
do
case $i in
    -g=*|--genome=*)
    GENOME="${i#*=}" #hg19
    shift # past argument=value
    ;;
    -s=*|--sampleFile=*)
    FILE_SAMPLE_LIST="${i#*=}" #test_sra_config.txt
    shift # past argument=value
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"  # /data/aryee/ellisen/tnbc/out_gencodeV24
    shift # past argument=value
    ;;
    -r=*|--strandedness=*)
    STRANDED="YES"
    STRANDEDNESS="${i#*=}"  # for hisat2  --rna-strandness FR or RF. Default unstranded
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done


if [[ ${GENOME} == "mm10" ]]; then
	GENOME_FASTA=/data/rivera/genomes/mm10/mm10.fa
	ANNOTATION_GTF=/data/rivera/genomes/mm10/gencode.vM6.basic.annotation.gtf
	if [[ ! -f  /data/rivera/sowmya/genomes/mm10.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes mm10  > /data/rivera/sowmya/genomes/mm10.chrom.sizes
	fi
elif [[ ${GENOME} == "hg19" ]]; then
	GENOME_FASTA=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
        ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf
	if [[ ! -f  /data/rivera/sowmya/genomes/hg19.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes hg19 > /data/rivera/sowmya/genomes/hg19.chrom.sizes
	fi
elif [[ ${GENOME} == "Zv9" ]]; then
	GENOME_FASTA=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
        ANNOTATION_GTF=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf 
	if [[ ! -f  /data/rivera/sowmya/genomes/Zv9.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes danRer7 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/Zv9.chrom.sizes
	fi
fi
echo $GENOME_FASTA $ANNOTATION_GTF

mkdir -p ${OUTPUT_DIR}/hisat2out
mkdir -p ${OUTPUT_DIR}/cufflinks_out



HISAT2_INDEX=/PHShome/si992/singlecell/txt/hg19/genome
HISAT2_SPLICE_SITES=$HOME/singlecell/txt/hisat2_splicesites_gencode_v24.txt 
ANNOTATION_GTF=/PHShome/si992/singlecell/txt/gencode.v24lift37.basic.annotation.gtf

echo "starting hisat lines"
while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`
	SRA=`echo ${line_all} | cut -d" " -f2`

        echo sample name is $SAMPLE_NAME

	hisat2_command="hisat2 -p 12 --dta --known-splicesite-infile ${HISAT2_SPLICE_SITES} -x ${HISAT2_INDEX} --sra-acc ${SRA} -S ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sam 2>${OUTPUT_DIR}/hisat2out/hisat2_log.${SAMPLE_NAME}.txt"

	echo """
        module load python/2.7.3
        module load RSeQC/2.4
        module load ucsc
	module load cufflinks/2.2.1

        which hisat2
        ${hisat2_command}
        samtools view -bS -o ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.raw.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sam
        samtools view -b -F 4 -o ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.bam  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.raw.bam
        samtools sort -T ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted -o ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.bam  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.bam
        java  -Xmx40g -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.bam O=${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.bam M=${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

        split_bam.py -i ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.bam -r ${HOME}/commonscripts/rnaseq/${GENOME}_rRNA.bed --out-prefix=${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed > ${OUTPUT_DIR}/hisat2out/rRNA_${SAMPLE_NAME}.txt
        mv ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.ex.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.bam
        samtools index ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.bam
        rm ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.in.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.junk.bam

        samtools flagstat ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.raw.bam >  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.flagstat
        totalRawReads=\`awk '{ if (NR == 1) print \$1 }'  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.flagstat\`
        mappedReads=\`awk '{ if (NR == 5) print \$1 }'  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.flagstat\`
        dedupedReads=\`samtools view -c  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.bam\`
        fraction_dups=\`grep \"Unknown \" ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.dups.txt | awk '{ print \$(NF-1)}'\`
        reads_after_rRNA_removal=\`samtools view -c  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.bam \`
        rm ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.raw.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.bam

        echo ${SAMPLE_NAME} \${totalRawReads} \${mappedReads} \${dedupedReads} \${fraction_dups} \${reads_after_rRNA_removal} > ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.mapping_metrics.txt 
	cufflinks -u -o ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME} -G ${ANNOTATION_GTF} ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.deduped.rRNA_removed.bam
	
	""" > ../bsubFiles/hisat2_alignment_${GENOME}_${SAMPLE_NAME}.bsub
done < ${FILE_SAMPLE_LIST}
