#GENOME=$1 #Zv9
#FILE_SAMPLE_LIST=$2 #test_Zv9_config.txt
#OUTPUT_DIR=$3 # /data/langenau/out
#COMBINE_FASTQ=$4 # true or false


if [[ $# -lt 3 ]]; then
    echo "Usage: `basename $0` -g=genome -s=sampleFile -o=outputDir --combineFastq(if lanes are to be combined)... ">&2
    exit 1
fi

COMBINE_FASTQ="NO"

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
    --combineFastq)
    COMBINE_FASTQ="YES"
    shift # past argument with no value
    ;;
    --help)
    echo "Usage: `basename $0` -g=genome -s=sampleFile -o=outputDir --combineFastq(if lanes are to be combined)... ">&2
    exit 0
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
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/hg19
	if [[ ! -f  /data/rivera/sowmya/genomes/hg19.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes hg19 > /data/rivera/sowmya/genomes/hg19.chrom.sizes
	fi
elif [[ ${GENOME} == "Zv9" ]]; then
	GENOME_FASTA=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
	BWA_INDEX=/PHShome/si992/commonscripts/bwa_index/Zv9
	if [[ ! -f  /data/rivera/sowmya/genomes/Zv9.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes danRer7 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/Zv9.chrom.sizes
	fi
fi
echo $GENOME_FASTA $BWA_INDEX $ANNOTATION_GTF

mkdir -p ${OUTPUT_DIR}/fastqc_out
mkdir -p ${OUTPUT_DIR}/bwa_out
mkdir -p ${OUTPUT_DIR}/bigwigs
mkdir -p ${OUTPUT_DIR}/combined_fastqs

while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`
	if [[ ${COMBINE_FASTQ} == "YES" ]]; then
	        FASTQ_PREFIX=`echo ${line_all} | cut -d" " -f2`
		R1_FASTQ=`\ls ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz`
        	R2_FASTQ=`\ls ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz`
		readLength=`zcat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz | head -2 | tail -1 | awk '{ print length($0)}'` # only checks read length of first read. Maybe a problem if reads have been trimmed or have variable lengths for other reasons
		echo """
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
		
		""" > ../bsubFiles/runFastqc_${SAMPLE_NAME}.bsub
	else
        	R1_FASTQ=`echo ${line_all} | cut -d" " -f2`
	        R2_FASTQ=`echo ${line_all} | cut -d" " -f3`
		ls ${R1_FASTQ} ${R2_FASTQ}
		readLength=`zcat ${R1_FASTQ} | head -2 | tail -1 | awk '{ print length($0)}'` # only checks read length of first read. Maybe a problem if reads have been trimmed or have variable lengths for other reasons.
		echo """
                        module load fastqc
                        fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out ${R1_FASTQ}
                        fastqc --noextract -o ${OUTPUT_DIR}/fastqc_out ${R2_FASTQ}
                """ > ../bsubFiles/runFastqc_${SAMPLE_NAME}.bsub
	fi

done < ${FILE_SAMPLE_LIST}


while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`

	if [[ ${COMBINE_FASTQ} == "YES" ]]; then
                FASTQ_PREFIX=`echo ${line_all} | cut -d" " -f2`
	fi

        R1_FASTQ=`echo ${line_all} | cut -d" " -f2`
        R2_FASTQ=`echo ${line_all} | cut -d" " -f3`

        echo sample name is $SAMPLE_NAME
        echo """
        module load bwa
        module load ucsc
	""" > ../bsubFiles/bwa_alignment_${GENOME}_${SAMPLE_NAME}.bsub

	if [[ ${COMBINE_FASTQ} == "YES" ]]; then
                FASTQ_PREFIX=`echo ${line_all} | cut -d" " -f2`
		echo """
	cat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
	cat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz
		""" >> ../bsubFiles/hisat2_alignment_${GENOME}_${SAMPLE_NAME}.bsub
		R1_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
		R2_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz
        fi

	echo """
	bwa aln -t 4 ${BWA_INDEX} ${R1_FASTQ} 2>${LOG_FILE_BWA} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai
	bwa aln -t 4 ${BWA_INDEX} ${R2_FASTQ} 2>${LOG_FILE_BWA} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai
	bwa sampe ${BWA_INDEX} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${R1_FASTQ} ${R2_FASTQ} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam
	samtools view -bS ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam
	samtools view -b -F2304 -q ${MIN_QUAL} ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
        samtools sort -T ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted -o ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.bam
        java  -Xmx40g -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.bam O=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam M=${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR

        samtools index ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam
        rm ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r1.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.r2.sai ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sam 

        samtools flagstat ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.raw.bam >  ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.flagstat
        totalRawReads=\`awk '{ if (NR == 1) print \$1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat\`
        mappedReads=\`awk '{ if (NR == 5) print \$1 }'  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.flagstat\`
        dedupedReads=\`samtools view -c  ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.sorted.deduped.bam\`
        fraction_dups=\`grep \"Unknown \" ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.dups.txt | awk '{ print \$(NF-1)}'\`
        rm ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.raw.bam ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sorted.bam

        echo ${SAMPLE_NAME} \${totalRawReads} \${mappedReads} \${dedupedReads} \${fraction_dups} > ${OUTPUT_DIR}/bwa_out/${SAMPLE_NAME}.mapping_metrics.txt """ >> ../bsubFiles/bwa_alignment_${GENOME}_${SAMPLE_NAME}.bsub
done < ${FILE_SAMPLE_LIST}
