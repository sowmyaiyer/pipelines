OUTPUT_DIR=$1 

cat ${OUTPUT_DIR}/STAR_out/*.mapping_metrics.txt > ${OUTPUT_DIR}/STAR_out/mapping_metrics_noheader.txt
echo "sample totalFragments uniquelyMappingReads multimappingReads fractionDup rRNA_removed_reads"> ${OUTPUT_DIR}/STAR_out/mapping_metrics_header
cat ${OUTPUT_DIR}/STAR_out/mapping_metrics_header ${OUTPUT_DIR}/STAR_out/mapping_metrics_noheader.txt | sed 's/ /\t/g' > ${OUTPUT_DIR}/STAR_out/mapping_metrics.txt

#Rscript plot_mapping_metrics.R ${OUTPUT_DIR}
