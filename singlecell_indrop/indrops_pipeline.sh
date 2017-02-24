#python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml build_index     \
#	--genome-fasta-gz /data/langenau/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz     \
#	--ensembl-gtf-gz /data/langenau/Danio_rerio.GRCz10.85.gtf.gz

# RSEM can not recognize transcript ENSDART00000000992|GTSF1!
# Gene Names with spaces in the GFT file were a problem with RSEM. So had to change (1 of many) to _1_of_many in GTF, gzip it and recreate index
#python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml build_index     \
#        --genome-fasta-gz /data/langenau/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz     \
#        --ensembl-gtf-gz /data/langenau/Danio_rerio.GRCz10.85.gtf.quoted_genenames.gz


# Round 2 of index building. Noticed that indrops.py was throwing out a lot of biotypes that were of interest to us. So changed that, reran indexing; reran counting
# STEP 0 - ONE TIME ONLY
python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml build_index     \
        --genome-fasta-gz /data/langenau/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz     \
        --ensembl-gtf-gz /data/langenau/Danio_rerio.GRCz10.85.gtf.quoted_genenames.gz

# STEP 1
for library in {"WT_3","PRKDC_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
do
        for worker_index in {0..39}
        do
		bsub -q big-multi python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml filter --libraries ${library} --total-workers 40 --worker-index ${worker_index}
	done
done

# STEP 2
python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml identify_abundant_barcodes
grep -c "\-" ./*/abundant_barcodes.pickle


# STEP 3
for library in {"WT_3","PRKDC_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
do
        for worker_index in {0..39}
        do
		bsub -q big-multi -o $HOME/RECYCLE_BIN/sort_worker_${worker_index}.log python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml sort --libraries ${library} --total-workers 40 --worker-index ${worker_index}
	done
done


# Had to use samtools/1.3.1 for samtools sort to work without changing code in indrops.py. Changed setting in /PHShome/si992/packages/indrops/test/indrops_v3.yaml
# module load samtools/1.3.1; which samtools;/apps/lib-osver/samtools/1.3.1/bin/. But this did not work!! 
# So I had to install samtools 1.3.1 in my home directory, change setting again in the yaml file to /PHShome/si992/packages/samtools-1.3.1/ and rerun
# STEP 4
for library in {"WT_3","PRKDC_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
do
	for worker_index in {0..39}
	do
		bsub -q big-multi -o $HOME/RECYCLE_BIN/quantify_worker_${worker_index}.log python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml quantify --min-reads 10000 --min-counts 1 --analysis-prefix round2.12_19_2016 --libraries ${library} --total-workers 40 --worker-index ${worker_index}
	done
done
# STEP 5
for library in {"WT_3","PRKDC_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
do
	bsub -q big-multi python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml quantify --analysis-prefix round2.12_19_2016  --total-workers 1 --libraries ${library}
done


# STEP 6 construct bigwigs
OUTDIR="/data/langenau/indrops_out/QT_inDrop_1205_and_06_2016/"
mkdir -p ${OUTDIR}/bigwigs
#for library in {"WT_3","PRKDC_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
for library in {"WT_3","IL2RGA_1","IL2RGA_2","PRKDC_and_IL2RGA_1","PRKDC_and_IL2RGA_2"}
do
	echo "start bw file generation"
	reads_total=`samtools view -c ${OUTDIR}/${library}/quant_dir/round2.12_19_2016.worker0_1.bam`
	scalingFactor=`echo ${reads_total} | awk -vmappedReads=${reads_total} 'BEGIN{print 1000000/mappedReads}'`
	echo ${reads_total} $scalingFactor
	bedtools genomecov -ibam ${OUTDIR}/${library}/quant_dir/round2.12_19_2016.worker0_1.bam -g  /data/rivera/sowmya/genomes/GRCz10.chrom.sizes -bga -split -scale ${scalingFactor} >  ${OUTDIR}/bigwigs/${library}.sorted.deduped.bg
	bedGraphToBigWig ${OUTDIR}/bigwigs/${library}.sorted.deduped.bg  /data/rivera/sowmya/genomes/GRCz10.chrom.sizes ${OUTDIR}/bigwigs/${library}.bw
	rm  ${OUTDIR}/bigwigs/${library}.sorted.deduped.bg
echo "DONE bw file generation for ${library}"
done
