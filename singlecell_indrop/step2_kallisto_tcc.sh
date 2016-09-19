#module load cufflinks
#gffread -w /data/langenau/Danio_rerio.Zv9.77.gtf.fa  -g /data/langenau/Danio_rerio.Zv9.77.dna.toplevel.fa /data/langenau/Danio_rerio.Zv9.77.gtf
#NOT USED kallisto index --make-unique -i ~/commonscripts/kallisto_index_Zv9 /data/langenau/Danio_rerio.Zv9.77.gtf.fa
#NOT USED kallisto index -i /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa.kallisto_index /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
#NOT USED kallisto index -i /data/langenau/Danio_rerio.GRCz10.cdna.all.fa.kallisto_index /data/langenau/Danio_rerio.GRCz10.cdna.all.fa

for INPUT_FASTQ_GZ in /data/langenau/singlecell_prkdc/fastq/Zebrafish*.processed.fastq.gz
do 
	FNAME=`basename ${INPUT_FASTQ_GZ} | sed 's/\.fastq\.gz//g'`
	echo """
	kallisto pseudo -i /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa.kallisto_index -o /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/${FNAME} --umi -b /data/langenau/singlecell_prkdc/demultiplexed/${FNAME}/barcodes.batch -t 4
	umis kallisto_table /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/${FNAME} /data/langenau/Danio_rerio.GRCz10.cdna.all.fa
	""" > ../bsubFiles/tcc_kallisto_${FNAME}.bsub
done 
