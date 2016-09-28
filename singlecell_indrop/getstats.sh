for sample in {"Zebrafish-WT1","Zebrafish-WT2","Zebrafish-Mut1_S1","Zebrafish-PRKDC2"}
do
	echo $sample
	if [[ -f  /data/langenau/singlecell_prkdc/out_singlecell/stats.${sample}.3.txt ]]; then
		rm  /data/langenau/singlecell_prkdc/out_singlecell/stats.${sample}.3.txt
	fi
	echo barcode reads mapped_reads total_umis detected_transcripts > /data/langenau/singlecell_prkdc/out_singlecell/stats.${sample}.3.txt
	while read bc 
	do
		echo $bc
		total_reads=`zcat /data/langenau/singlecell_prkdc/demultiplexed/${sample}.processed/$bc.fastq.gz | wc -l`
		total_umis=`awk -F"," 'BEGIN{sum=0}{sum=sum+$2}END{print sum}' /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.${bc}.genename.umi_counts`
		detected_transcripts=`awk -F"," 'BEGIN{sum=0}{if ($2 > 0) {sum=sum+1}}END{print sum}' /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.${bc}.genename.umi_counts`
		echo $bc $((total_reads/4)) `cut -f1 /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.${bc}.OUT.sam | sort | uniq | wc -l` ${total_umis} ${detected_transcripts} >> /data/langenau/singlecell_prkdc/out_singlecell/stats.${sample}.3.txt
	done < /data/langenau/singlecell_prkdc/demultiplexed/${sample}.processed.barcodes_above_threshold.txt
done
