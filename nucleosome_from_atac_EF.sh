# Remove PCR dups. 
java -jar /apps/source/picardtools/picard-tools-1.58/MarkDuplicates.jar I=/PHShome/gz976/rivera/ewings/ATACseqProcessingNew/ATAC.combi.MSC.EF.nodups.bam O=/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.bam REMOVE_DUPLICATES=true M=/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.txt
samtools index /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.bam
samtools view -b /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.properchrs.bam

# Get properly paired reads and separate reads by fragment size - thresholds from Greenleaf paper
samtools view -h /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.properchrs.bam | awk  '{
                 if ($2 == 83 || $2 == 163 || $2 == 99 || $2 == 147)
                 {
                         isize = sqrt($9*$9)
                         if (isize < 100)
                                 print $0 >> "/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.tf.sam"
                         else if (isize >= 180 && isize <= 247)
                                 print $0 >> "/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.mononuc.sam"
                         else if (isize >= 315 && isize <= 473)
                                 print $0 >> "/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.dinuc.sam"
                         else if (isize >= 558 && isize <= 615)
                                 print $0 >> "/PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.trinuc.sam"
                 }
 }'
samtools view -h /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.properchrs.bam | grep -P "^@" > /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.header_only.sam
cat /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.header_only.sam /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.tf.sam | samtools view -b - | samtools sort -n -T ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf -o ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.bam
cat /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.header_only.sam /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.mononuc.sam | samtools view -b - | samtools sort -n -T ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mono -o ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.bam
cat /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.header_only.sam /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.dinuc.sam | samtools view -b - | samtools sort -n -T ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.di -o ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.bam
cat /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.header_only.sam /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.trinuc.sam | samtools view -b - | samtools sort -n -T ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tri -o ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.bam

bedtools bamtobed -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bed
bedtools bamtobed -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.bed
bedtools bamtobed -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.fragments.bed
bedtools bamtobed -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.bam -bedpe | cut -f1,2,6 | sort -k1,1 -k2,2n > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.fragments.bed

# Split dinucleotide fragments into two mono fragments 
awk '{ 
	len=$3-$2; 
	half=len/2; 
	printf("%s\t%d\t%d\n%s\t%d\t%d\n", $1,$2,$2+half,$1,$2+half+1,$3)
}' ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.fragments.bed > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.fragments.mono.bed

awk '{
        len=$3-$2;
        third=len/3;
	two_thirds=2*len/3
        printf("%s\t%d\t%d\n%s\t%d\t%d\n%s\t%d\t%d\n", $1,$2,$2+third,$1,$2+third+1,$2+two_thirds,$1,$2+two_thirds+1,$3)
}' ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.fragments.bed > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.fragments.mono.bed


cat ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.bed ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.dinuc.fragments.mono.bed ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.trinuc.fragments.mono.bed | sort -k1,1 -k2,2n > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.all.bed

# See if we can get a signal from tf and mono fragments. Also get signal from all atac fragments
read_depth=`samtools view -c /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.properchrs.bam`
scale_factor=`echo $read_depth | awk -vmappedReads=$read_depth 'BEGIN{print 10000000/mappedReads}'`
bedtools genomecov -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.all.bed -g ../txt/hg19.nochr.len -bga -scale ${scale_factor}  > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.all.bedGraph
bedGraphToBigWig ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.all.bedGraph ../txt/hg19.nochr.len ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.all.bw

bedtools genomecov -i ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bed -g ../txt/hg19.nochr.len -bga -scale ${scale_factor} > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bedGraph
bedGraphToBigWig ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bedGraph ../txt/hg19.nochr.len ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bw

bedtools genomecov -ibam /PHShome/si992/baf_ews/txt/ATAC.combi.MSC.EF.nodups.nodups.properchrs.bam -g ../txt/hg19.nochr.len -bga -scale ${scale_factor} > ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.bedGraph
bedGraphToBigWig ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.bedGraph ../txt/hg19.nochr.len ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.bw


/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.tf.fragments.bw output.tf.EF.txt
/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.mononuc.fragments.all.bw output.mono.EF.txt
/apps/lab/aryee/bwtool/bin/bwtool matrix 3000:3000 -tiled-averages=10 /PHShome/gz976/rivera/ewings/ConsensusRepeats.NOchr.txt ../txt/ATAC.combi.MSC.EF.nodups.nodups.sorted_by_name.bw output.allfragments.EF.txt
