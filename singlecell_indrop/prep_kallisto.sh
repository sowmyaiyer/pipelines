# Merged fa from /data/langenau/Danio_rerio.GRCz10.cdna.all.fa from ensembl ftp site ftp ftp.ensembl.org /pub/release-85/fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz and 
# /data/langenau/Danio_rerio.GRCz10.ncrna.fa from /pub/release-85/fasta/danio_rerio/ncrna/Danio_rerio.GRCz10.ncrna.fa.gz
# Had to do this because they had moved a T-cell receptor gene from protein-coding to lincRNA status, so it was in a different fa file

cat /data/langenau/Danio_rerio.GRCz10.cdna.all.fa /data/langenau/Danio_rerio.GRCz10.ncrna.fa > /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa
kallisto index -i /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa.kallisto_index /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa

grep ">ENS" /data/langenau/Danio_rerio.GRCz10.cdna_and_ncrna.fa | awk '{ 
							if ($2 == "cdna:novel") 
								print $1"\t"$4"\t"$4; 
							else print $1"\t"$4"\t"$7 
						}' | sed 's/^>//g' | sed 's/gene://g' | sed 's/gene_symbol://g' > ~/langenau/txt/ensemblTranscript_conversion_table_from_fa.txt

cut -f1,2 ~/langenau/txt/ensemblTranscript_conversion_table_from_fa.txt > ~/langenau/txt/ensemblTranscript_to_ensemblGene.from_fa.txt
cut -f1,3 ~/langenau/txt/ensemblTranscript_conversion_table_from_fa.txt > ~/langenau/txt/ensemblTranscript_to_gene_symbol.from_fa.txt
