for sample in "Zebrafish-WT1" "Zebrafish-WT2" "Zebrafish-PRKDC2" "Zebrafish-Mut1_S1"
do
	echo ${sample}
	if [[ -f ${sample}.processed.sampleFile ]]; then
		rm ${sample}.processed.sampleFile
	fi
	for f in /data/langenau/singlecell_prkdc/demultiplexed/${sample}.processed/*.fastq.gz
	do
		sample_name=`basename $f | sed 's/\.fastq\.gz//g'`
		echo ${sample}"."$sample_name $f >> ${sample}.processed.sampleFile
	done
done
