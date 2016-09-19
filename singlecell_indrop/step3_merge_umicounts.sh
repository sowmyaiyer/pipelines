# checking order of genes to paste. Take a random cell and compare genes from all other cells to it.
#cut -d"," -f1  /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/Zebrafish-Mut1_S1.AGCCAAGAT-GGGCCAAT.genename.umi_counts > ~/test1
#for sample in {"Zebrafish-WT1","Zebrafish-WT2","Zebrafish-Mut1_S1","Zebrafish-PRKDC2"}
#do
# 	for counts_file in /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.*.genename.umi_counts
#	do
#		echo $counts_file
#		cut -d"," -f1 $counts_file > ~/genenames.`basename $counts_file`.txt
#		diff ~/test1  ~/genenames.`basename $counts_file`.txt
#	done
#done
#cut -d"," -f1  /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/Zebrafish-Mut1_S1.AGCCAAGAT-GGGCCAAT.ensemblgenename.umi_counts > ~/test1.ensembl
#for sample in {"Zebrafish-WT1","Zebrafish-WT2","Zebrafish-Mut1_S1","Zebrafish-PRKDC2"}
#do
#        for counts_file in /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.*.ensemblgenename.umi_counts
#        do
#                echo $counts_file
#                cut -d"," -f1 $counts_file > ~/ensemblgenenames.`basename $counts_file`.txt
#                diff ~/test1.ensembl  ~/ensemblgenenames.`basename $counts_file`.txt
#        done
#done
# END checking order of genes to paste.

for sample in {"Zebrafish-WT1","Zebrafish-WT2","Zebrafish-Mut1_S1","Zebrafish-PRKDC2"}
do
	echo $sample
        paste /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.*.genename.umi_counts | awk -F"\t" '{
		split($1,arr1,",");
		printf("%s\t",arr1[1]);
		for (i = 2; i <= (NF-2); i++ )
		{	
			split($i,arr,",")
			printf("%s\t",arr[2])
		}
		split($NF, arr2, ",")
		printf("%s\n", arr2[2])
	}' > /data/langenau/singlecell_prkdc/out_singlecell/${sample}.all.genename.umi_counts
done

for sample in {"Zebrafish-WT1","Zebrafish-WT2","Zebrafish-Mut1_S1","Zebrafish-PRKDC2"}
do
        echo $sample
        paste /data/langenau/singlecell_prkdc/out_singlecell/out_kallisto/${sample}.*.ensemblgenename.umi_counts | awk -F"\t" '{
                split($1,arr1,",");
                printf("%s\t",arr1[1]);
                for (i = 2; i <= (NF-2); i++ )
                {
                        split($i,arr,",")
                        printf("%s\t",arr[2])
                }
                split($NF, arr2, ",")
                printf("%s\n", arr2[2])
        }' > /data/langenau/singlecell_prkdc/out_singlecell/${sample}.all.ensemblgenename.umi_counts
done
