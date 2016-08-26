awk '{ printf(">barcode1_%s_%s\n%s\n",NR,$0,$0) }' /data/langenau/singlecell_prkdc/gel_barcode1_list_revComp.txt > /data/langenau/singlecell_prkdc/gel_barcode1_list_revComp.fa
awk '{ printf(">barcode2_%s_%s\n%s\n",NR,$0,$0) }' /data/langenau/singlecell_prkdc/gel_barcode2_list_revComp.txt > /data/langenau/singlecell_prkdc/gel_barcode2_list_revComp.fa

grep "^@NB" Zebrafish-WT1.parsed.fastq |  cut -d":" -f8 | sed 's/CELL_//g' | cut -d"-" -f1  > debug_Zebrafish-WT1.parsed_bc1.txt

module load bowtie
bowtie-build /data/langenau/singlecell_prkdc/gel_barcode1_list_revComp.fa gel_barcode1_list_revComp.fa
bowtie-build /data/langenau/singlecell_prkdc/gel_barcode2_list_revComp.fa gel_barcode2_list_revComp.fa
#bowtie -v 3 --all --norc -S /PHShome/si992/commonscripts/singlecell_indrop/gel_barcode1_list_revComp.fa -r /data/langenau/singlecell_prkdc/fastq/debug_Zebrafish-WT1.parsed_bc1.txt 
bowtie -v 3 --all --norc -S /PHShome/si992/commonscripts/singlecell_indrop/gel_barcode1_list_revComp.fa -f /data/langenau/singlecell_prkdc/fastq/debug_Zebrafish-WT1.parsed_bc1.fa > /data/langenau/singlecell_prkdc/fastq/debug_Zebrafish-WT1.parsed_bc1.out.sam
bowtie -v 3 --all --norc -S /PHShome/si992/commonscripts/singlecell_indrop/gel_barcode1_list_revComp.fa -f /data/langenau/singlecell_prkdc/fastq/debug_Zebrafish-WT1.parsed_bc1.fa > /data/langenau/singlecell_prkdc/fastq/debug_Zebrafish-WT1.parsed_bc1.out.sam

# Get edit distance between barcodes by mapping barcodes to themselves
bowtie -v 3 --all --norc -S /PHShome/si992/commonscripts/singlecell_indrop/gel_barcode1_list_revComp.fa -f ../gel_barcode1_list_revComp.fa  | awk '{ if ($1 != $3 && $4 == 1) print }'

# Got invalid barcodes list (/data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq.txt) from ~/langenau/scripts/indrop_reverse_barcode.rmd
echo """
module load bowtie
bowtie -v 3 --all --norc -S /PHShome/si992/commonscripts/singlecell_indrop/gel_barcode1_list_revComp.fa -r /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq.txt > /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq_mapped_to_bc1.sam
samtools view -b /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq_mapped_to_bc1.sam > /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq_mapped_to_bc1.bam
rm /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq_mapped_to_bc1.sam
""" > ../bsubFiles/map_invalid_barcode1.bsub

samtools view  -F4 invalid_bc1s_in_fastq_mapped_to_bc1.bam | awk '{ split($3,arr,"_"); if (($4 == 1) && length(arr[3]) == length($10)) print }' > invalid_bc1s_in_fastq_mapped_to_bc1.mismatches.txt


sort /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq.txt | uniq -c | awk '{ print $2"\t"$1}' | sort -k2,2nr > /data/langenau/singlecell_prkdc/invalid_bc1s_in_fastq.dist.txt



# Get filtered out reads and see how their barcodes map to known barcodes


