GENOME_FASTA=/pub/genome_references/UCSC/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
STAR_INDEX_DIR=/data/rivera/sowmya/genomes/STAR_indices/STAR_index_hg19
ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.gtf
echo """
module load aryee/star-2.4.0h
STAR --runMode genomeGenerate --genomeDir $STAR_INDEX_DIR --genomeFastaFiles $GENOME_FASTA --sjdbGTFfile $ANNOTATION_GTF  --sjdbOverhang 99 --runThreadN 12
""" > ../bsubFiles/build_hg19_STAR_index.bsub
