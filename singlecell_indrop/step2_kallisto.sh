#module load cufflinks
#NOT USED gffread -w /data/langenau/Danio_rerio.Zv9.77.gtf.fa  -g /data/langenau/Danio_rerio.Zv9.77.dna.toplevel.fa /data/langenau/Danio_rerio.Zv9.77.gtf
#NOT USED kallisto index --make-unique -i ~/commonscripts/kallisto_index_Zv9 /data/langenau/Danio_rerio.Zv9.77.gtf.fa
#NOT USED kallisto index -i /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa.kallisto_index /data/langenau/Danio_rerio.Zv9.rel79.cdna.all.fa
#NOT USED kallisto index -i /data/langenau/Danio_rerio.GRCz10.cdna.all.fa.kallisto_index /data/langenau/Danio_rerio.GRCz10.cdna.all.fa
# SEE ./buildGeneNameTable.sh

while read lineall
do 
	sample_name=`echo $lineall | awk '{ print $1}'`
	fastq=`echo $lineall | awk '{ print $2}'`
	echo """
	kallisto quant --pseudobam --single -l 200 -s 20 -i /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa.kallisto_index -o /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name} ${fastq} > /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name}.OUT.sam
	umis tagcount --genemap $HOME/langenau/txt/ensemblTranscript_to_gene_symbol.from_fa.txt /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name}.OUT.sam /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name}.genename.umi_counts
	umis tagcount --genemap $HOME/langenau/txt/ensemblTranscript_to_ensemblGene.from_fa.txt /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name}.OUT.sam /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample_name}.ensemblgenename.umi_counts
	""" > ../bsubFiles/kallisto_${sample_name}.bsub
done < Zebrafish-all.sampleFile
