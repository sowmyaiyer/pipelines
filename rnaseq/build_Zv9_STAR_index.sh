GENOME_FASTA=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
ANNOTATION_GTF=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf
STAR_INDEX_DIR=/PHShome/si992/commonscripts/rnaseq/STAR_indices/STAR_index_Zv9
echo """
module load aryee/star-2.4.0h
STAR --runMode genomeGenerate --genomeDir $STAR_INDEX_DIR --genomeFastaFiles $GENOME_FASTA --sjdbGTFfile $ANNOTATION_GTF --sjdbOverhang 99 --runThreadN 12
""" > ../bsubFiles/build_Zv9_STAR_index.bsub
