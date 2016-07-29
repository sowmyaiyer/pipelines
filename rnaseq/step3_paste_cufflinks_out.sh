OUTPUT_DIR=$1 # /data/rivera/lombard/out
SAMPLE_FILE=$2 # $HOME/commonscripts/rnaseq/test_lombard_config.txt
while read line_all
do
	SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`
	sed '1d' ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME}/genes.fpkm_tracking | awk -F"\t" -vsamplename=${SAMPLE_NAME} 'BEGIN{print "gene\tlocus\t"samplename} {print $1"\t"$7"\t"$10}' > ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME}_fpkm_with_header.txt
done < ${SAMPLE_FILE}


paste ${OUTPUT_DIR}/cufflinks_out/*_fpkm_with_header.txt > ${OUTPUT_DIR}/cufflinks_out/fpkm_all_genes.txt

awk '{
        printf("%s\t%s\t%s\t",$1,$1,$2);
        for (i = 1; i < NF; i += 3)
        {
                if ($1 == $i && $2 == $(i+1))
                        printf ("%s\t",$(i+2))
        }
       printf("\n")
}' ${OUTPUT_DIR}/cufflinks_out/fpkm_all_genes.txt > ${OUTPUT_DIR}/cufflinks_out/fpkm_all_genes.final.txt


while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`
        sed '1d' ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME}/isoforms.fpkm_tracking | awk -F"\t" -vsamplename=${SAMPLE_NAME} 'BEGIN{print "gene\tlocus\t"samplename} {print $1"\t"$7"\t"$10}' > ${OUTPUT_DIR}/cufflinks_out/${SAMPLE_NAME}_fpkm_with_header.txt
done < ${SAMPLE_FILE}


paste ${OUTPUT_DIR}/cufflinks_out/*_fpkm_with_header.txt > ${OUTPUT_DIR}/cufflinks_out/fpkm_all_isoforms.txt

awk '{
        printf("%s\t%s\t%s\t",$1,$1,$2);
        for (i = 1; i < NF; i += 3)
        {
                if ($1 == $i && $2 == $(i+1))
                        printf ("%s\t",$(i+2))
        }
       printf("\n")
}'  ${OUTPUT_DIR}/cufflinks_out/fpkm_all_isoforms.txt > ${OUTPUT_DIR}/cufflinks_out/fpkm_all_isoforms.final.txt
