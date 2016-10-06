INPUT_BAM_FILE=$1 # /PHShome/gz976/rivera/ewings/ATACseqProcessingNew/ATAC.combi.MSC.EF.nodups.bam
SAMPLE_NAME=`basename $INPUT_BAM_FILE | sed 's/.bam//g'`
OUTPUT_DIR=$2 # /PHShome/si992/baf_ews/txt
OUTPUT_PREFIX=${OUTPUT_DIR}/${SAMPLE_NAME}


# Separate reads by fragment size - thresholds from Greenleaf paper
rm ${OUTPUT_PREFIX}*sam
samtools view -h ${INPUT_BAM_FILE} | awk  -vOUTPUT_PREFIX=${OUTPUT_PREFIX} '{
                 if ($2 == 83 || $2 == 163 || $2 == 99 || $2 == 147)
                 {
                         isize = sqrt($9*$9)
                         if (isize < 100)
                                 print $0 >> OUTPUT_PREFIX".nodups.nodups.tf.sam"
                         else if (isize >= 180 && isize <= 247)
                                 print $0 >> OUTPUT_PREFIX".nodups.nodups.mononuc.sam"
                         else if (isize >= 315 && isize <= 473)
                                 print $0 >> OUTPUT_PREFIX".nodups.nodups.dinuc.sam"
                         else if (isize >= 558 && isize <= 615)
                                 print $0 >> OUTPUT_PREFIX".nodups.nodups.trinuc.sam"
                 }
 }'
samtools view -H  ${INPUT_BAM_FILE} > ${OUTPUT_PREFIX}.nodups.nodups.header_only.sam
cat ${OUTPUT_PREFIX}.nodups.nodups.header_only.sam ${OUTPUT_PREFIX}.nodups.nodups.tf.sam | samtools view -b - | samtools sort -n -T ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf -o ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.bam
cat ${OUTPUT_PREFIX}.nodups.nodups.header_only.sam ${OUTPUT_PREFIX}.nodups.nodups.mononuc.sam | samtools view -b - | samtools sort -n -T ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mono -o ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.bam
cat ${OUTPUT_PREFIX}.nodups.nodups.header_only.sam ${OUTPUT_PREFIX}.nodups.nodups.dinuc.sam | samtools view -b - | samtools sort -n -T ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.di -o ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.bam
cat ${OUTPUT_PREFIX}.nodups.nodups.header_only.sam ${OUTPUT_PREFIX}.nodups.nodups.trinuc.sam | samtools view -b - | samtools sort -n -T ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tri -o ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.bam

bedtools bamtobed -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bed
bedtools bamtobed -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.bed
bedtools bamtobed -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.fragments.bed
bedtools bamtobed -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.fragments.bed

# Split dinucleotide fragments into two mono fragments 
awk '{ 
	len=$3-$2; 
	half=len/2; 
	printf("%s\t%d\t%d\n%s\t%d\t%d\n", $1,$2,$2+half,$1,$2+half+1,$3)
}' ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.fragments.bed > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.fragments.mono.bed

awk '{
        len=$3-$2;
        third=len/3;
	two_thirds=2*len/3
        printf("%s\t%d\t%d\n%s\t%d\t%d\n%s\t%d\t%d\n", $1,$2,$2+third,$1,$2+third+1,$2+two_thirds,$1,$2+two_thirds+1,$3)
}' ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.fragments.bed > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.fragments.mono.bed


cat ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.bed ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.dinuc.fragments.mono.bed ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.trinuc.fragments.mono.bed | sort -k1,1 -k2,2n > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.all.bed

# See if we can get a signal from tf and mono fragments. Also get signal from all atac fragments
read_depth=`samtools view -c ${INPUT_BAM_FILE}`
scale_factor=`echo $read_depth | awk -vmappedReads=$read_depth 'BEGIN{print 10000000/mappedReads}'`
bedtools genomecov -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.all.bed -g $HOME/commonscripts/hg19.nochr.len -bga -scale ${scale_factor}  > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.all.bedGraph
bedGraphToBigWig ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.all.bedGraph $HOME/commonscripts/hg19.nochr.len ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.all.bw

bedtools genomecov -i ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bed -g $HOME/commonscripts/hg19.nochr.len -bga -scale ${scale_factor} > ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bedGraph
bedGraphToBigWig ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bedGraph $HOME/commonscripts/hg19.nochr.len ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bw


#/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.tf.fragments.bw output.tf.EF.txt
#/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.mononuc.fragments.all.bw output.mono.EF.txt
#/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ${OUTPUT_PREFIX}.nodups.nodups.sorted_by_name.bw output.allfragments.EF.txt
