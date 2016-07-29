organism=$1
if [[ ${organism} == "hg19" ]]; then
	GENOME=/pub/home/aryee/genomes/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
	ANNOTATION_GTF=/data/rivera/genomes/UCSC_refFLAT.07_08_2016/ucsc_refFlat.07_08_2016.NOCHR.gtf
elif [[ ${organism} == "mm10" ]]; then
	GENOME=/data/rivera/genomes/mm10/mm10.fa
	ANNOTATION_GTF=/data/rivera/genomes/mm10/gencode.vM6.basic.annotation.gtf
elif [[ ${organism} == "Zv9" ]]; then
	GENOME=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.dna.toplevel.fa
	ANNOTATION_GTF=/pub/genome_references/Zv9/Danio_rerio.Zv9.69.gtf
fi
echo ${GENOME} ${ANNOTATION_GTF}
echo """
source /apps/lab/aryee/pyenv/versions/venv-2.7.10/bin/activate
hisat2_extract_splice_sites.py ${ANNOTATION_GTF} > /PHShome/si992/commonscripts/rnaseq/hisat2_splicesites/${organism}_hisat2_splicesites.txt
hisat2-build ${GENOME} /PHShome/si992/commonscripts/rnaseq/hisat2_indices/${organism}
""" > ../bsubFiles/build_hisat2_index.${organism}.bsub
