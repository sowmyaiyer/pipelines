#python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml build_index     \
#	--genome-fasta-gz /data/langenau/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz     \
#	--ensembl-gtf-gz /data/langenau/Danio_rerio.GRCz10.85.gtf.gz

# RSEM can not recognize transcript ENSDART00000000992|GTSF1!
# Gene Names with spaces in the GFT file were a problem with RSEM. So had to change (1 of many) to _1_of_many in GTF, gzip it and recreate index
python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml build_index     \
        --genome-fasta-gz /data/langenau/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz     \
        --ensembl-gtf-gz /data/langenau/Danio_rerio.GRCz10.85.gtf.quoted_genenames.gz


for worker_index in {0..39}
do
bsub -q big-multi python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml filter --total-workers 40 --worker-index ${worker_index}
done

python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml identify_abundant_barcodes

for worker_index in {0..39}
do
bsub -q big-multi python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml sort --total-workers 40 --worker-index ${worker_index}
done


# Had to use samtools/1.3.1 for samtools sort to work without changing code in indrops.py. Changed setting in /PHShome/si992/packages/indrops/test/indrops_v3.yaml
# module load samtools/1.3.1; which samtools;/apps/lib-osver/samtools/1.3.1/bin/. But this did not work!! 
# So I had to install samtools 1.3.1 in my home directory, change setting again in the yaml file to /PHShome/si992/packages/samtools-1.3.1/ and rerun
for worker_index in {0..39}
do
bsub -q big-multi python $HOME/packages/indrops/indrops.py /PHShome/si992/packages/indrops/test/indrops_v3.yaml quantify --min-reads 10000 --min-counts 1 --analysis-prefix round1.12_19_2016 --total-workers 40 --worker-index ${worker_index}
done
