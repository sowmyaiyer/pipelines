#GENOME=$1 #Zv9
#FILE_SAMPLE_LIST=$2 #test_Zv9_config.txt
#OUTPUT_DIR=$3 # /data/langenau/out
#COMBINE_FASTQ=$4 # true or false


if [[ $# -lt 3 ]]; then
    echo "Usage: `basename $0` -g=genome -s=sampleFile -o=outputDir -r=(if stranded library, forward or reverse; if unstranded, do not provide) --combineFastq(if lanes are to be combined)... ">&2
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
	GENOME_FASTA=/data/langenau/Danio_rerio.Zv9.77.dna.toplevel.fa #/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
	STAR_INDEX_DIR=/data/langenau/STAR_index_Zv9.77 #/data/langenau/STAR_index_Zv9
        ANNOTATION_GTF=/data/langenau/Danio_rerio.Zv9.77.gtf #/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf 
	if [[ ! -f  /data/rivera/sowmya/genomes/Zv9.chrom.sizes ]]; then
		/data/aryee/pub/genomes/fetchChromSizes danRer7 | sed 's/^chr//g' > /data/rivera/sowmya/genomes/Zv9.chrom.sizes
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

if [[ ${STRANDEDNESS} == "reverse" ]]; then
        CUFFLINKS_STRANDEDNESS="fr-firststrand"
	HISAT2_STRANDEDNESS="RF"
elif [[ ${STRANDEDNESS} == "forward" ]]; then
        CUFFLINKS_STRANDEDNESS="fr-secondstrand"
	HISAT2_STRANDEDNESS="FR"
else
        CUFFLINKS_STRANDEDNESS="fr-unstranded"
	HISAT2_STRANDEDNESS="unstranded"
fi

mkdir -p ${OUTPUT_DIR}/fastqc_out
mkdir -p ${OUTPUT_DIR}/STAR_out
mkdir -p ${OUTPUT_DIR}/hisat2out
mkdir -p ${OUTPUT_DIR}/cufflinks_out
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

        echo sample name is ${SAMPLE_NAME}
	echo fastq_1 is ${R1_FASTQ}
	echo fastq_2 is ${R2_FASTQ}

	if [[ ${STRANDED} == "NO" ]]; then
		added_star_option=" --outSAMstrandField intronMotif"
	else
		added_star_option=""	
	fi

	sjdbOverhang=$((readLength -1))
        star_fastq_R1=`echo $R1_FASTQ | sed 's/ /,/g'`
        star_fastq_R2=`echo $R2_FASTQ | sed 's/ /,/g'`

	echo """
	module load aryee/star-2.4.0h
	module load samtools/1.1

	module load python/2.7.3
        module load RSeQC/2.4
        module load ucsc
	module load cufflinks/2.2.1

	star_fastq_R1=`echo $R1_FASTQ | sed 's/ /,/g'`
	star_fastq_R2=`echo $R2_FASTQ | sed 's/ /,/g'`	
	STAR ${added_star_option} --genomeDir $STAR_INDEX_DIR --genomeFastaFiles $GENOME_FASTA --sjdbGTFfile $ANNOTATION_GTF --sjdbOverhang ${sjdbOverhang} --runThreadN 1 --outSAMtype BAM Unsorted --outFileNamePrefix  ${OUTPUT_DIR}/STAR_out/$SAMPLE_NAME --readFilesCommand zcat --readFilesIn $star_fastq_R1 $star_fastq_R2
	
	# Sorting in separate step because STAR sometimes crashes while sorting by coordinate
	
	samtools sort -T ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out -o ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.out.bam
	
	# Remove duplicates
	java  -Xmx40g -jar /apps/source/picard-tools-1.95/MarkDuplicates.jar I=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam O=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.bam M=${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}.dups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR
	
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

        scalingFactor=\`echo \${reads_after_rRNA_removal} | awk -vmappedReads=\${reads_after_rRNA_removal} 'BEGIN{print 1000000/mappedReads}'\`
        bedtools genomecov -ibam ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam -g  /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes -bga -split -scale ${scalingFactor} >  ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg
        bedGraphToBigWig ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg   /data/rivera/sowmya/genomes/${GENOME}.chrom.sizes ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bw
	rm ${OUTPUT_DIR}/bigwigs/${SAMPLE_NAME}.sorted.deduped.bg

	# Get FPKMs for isoforms
	cufflinks --quiet --library-type ${CUFFLINKS_STRANDEDNESS} -o ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME} -G ${ANNOTATION_GTF} ${OUTPUT_DIR}/STAR_out/${SAMPLE_NAME}Aligned.sortedByCoord.out.deduped.rRNA_removed.bam
	
	"""  > ../bsubFiles/star_alignment_${GENOME}_${SAMPLE_NAME}.bsub
done < ${FILE_SAMPLE_LIST}



HISAT2_INDEX=/PHShome/si992/commonscripts/rnaseq/hisat2_indices/${GENOME}
HISAT2_SPLICE_SITES=/PHShome/si992/commonscripts/rnaseq/hisat2_splicesites/${GENOME}_hisat2_splicesites.txt
echo "starting hisat lines"
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
        module load python/2.7.3
        module load RSeQC/2.4
        module load ucsc

        which hisat2
	""" > ../bsubFiles/hisat2_alignment_${GENOME}_${SAMPLE_NAME}.bsub

	if [[ ${COMBINE_FASTQ} == "YES" ]]; then
                FASTQ_PREFIX=`echo ${line_all} | cut -d" " -f2`
		echo """
	cat ${FASTQ_PREFIX}*_L00*_R1_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
	cat ${FASTQ_PREFIX}*_L00*_R2_001.fastq.gz > ${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz
		""" >> ../bsubFiles/hisat2_alignment_${GENOME}_${SAMPLE_NAME}.bsub
		R1_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R1.fastq.gz
		R2_FASTQ=${OUTPUT_DIR}/combined_fastqs/${SAMPLE_NAME}_R2.fastq.gz
        fi
	if [[ ${STRANDED} == "YES" ]]; then
	hisat2_command="hisat2 -p 12 --rna-strandness ${HISAT2_STRANDEDNESS} --dta --known-splicesite-infile ${HISAT2_SPLICE_SITES} -x ${HISAT2_INDEX} -1 ${R1_FASTQ} -2 ${R2_FASTQ} -S ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sam 2>${OUTPUT_DIR}/hisat2out/hisat2_log.${SAMPLE_NAME}.txt"
	else
	hisat2_command="hisat2 -p 12 --dta --known-splicesite-infile ${HISAT2_SPLICE_SITES} -x ${HISAT2_INDEX} -1 ${R1_FASTQ} -2 ${R2_FASTQ} -S ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.sam 2>${OUTPUT_DIR}/hisat2out/hisat2_log.${SAMPLE_NAME}.txt"
	fi

	echo """
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
        echo ${SAMPLE_NAME} \${totalRawReads} \${mappedReads} \${dedupedReads} \${fraction_dups} \${reads_after_rRNA_removal} > ${OUTPUT_DIR}/hisat2out/${SAMPLE_NAME}.mapping_metrics.txt """ >> ../bsubFiles/hisat2_alignment_${GENOME}_${SAMPLE_NAME}.bsub
done < ${FILE_SAMPLE_LIST}
