./combineLanes.sh

./step1_process_fastq.sh /data/langenau/singlecell_prkdc/sample_file.txt /data/langenau/singlecell_prkdc/fastq/


./test_demultiplex.java.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT1.processed.fastq.gz /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT1.cb_histogram.txt /data/langenau/singlecell_prkdc/demultiplexed/
./test_demultiplex.java.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT2.processed.fastq.gz /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT2.cb_histogram.txt /data/langenau/singlecell_prkdc/demultiplexed/
./test_demultiplex.java.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-PRKDC2.processed.fastq.gz /data/langenau/singlecell_prkdc/fastq/Zebrafish-PRKDC2.cb_histogram.txt /data/langenau/singlecell_prkdc/demultiplexed/
./test_demultiplex.java.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-Mut1_S1.processed.fastq.gz /data/langenau/singlecell_prkdc/fastq/Zebrafish-Mut1_S1.cb_histogram.txt /data/langenau/singlecell_prkdc/demultiplexed/



bsub -q big-multi ./demultiplex_step2.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT1.processed.fastq.gz
bsub -q big-multi ./demultiplex_step2.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT2.processed.fastq.gz
bsub -q big-multi ./demultiplex_step2.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-Mut1_S1.processed.fastq.gz
bsub -q big-multi ./demultiplex_step2.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-PRKDC2.processed.fastq.gz


./check_demultiplex.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT1.processed.fastq.gz
./check_demultiplex.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-WT2.processed.fastq.gz
./check_demultiplex.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-Mut1_S1.processed.fastq.gz
./check_demultiplex.sh /data/langenau/singlecell_prkdc/fastq/Zebrafish-PRKDC2.processed.fastq.gz

./buildSampleFile.sh

if [ -f Zebrafish-all.sampleFile ]; then
	rm Zebrafish-all.sampleFile
fi
cat *sampleFile > Zebrafish-all.sampleFile
./prep_kallisto.sh 

# For kallisto + umi counting
./step2_kallisto.sh
./step3_merge_umicounts.sh




# For kallisto with equivalence classes
./build_barcodes_batch_file.sh
./step2_kallisto_tcc.sh

umis kallisto_table /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/Zebrafish-Mut1_S1.processed /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
umis kallisto_table /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/Zebrafish-WT1.processed /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
umis kallisto_table /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/Zebrafish-WT2.processed /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
umis kallisto_table /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto_tcc/Zebrafish-PRKDC2.processed /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
