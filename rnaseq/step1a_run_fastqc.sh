FILE_SAMPLE_LIST=test_mm10_config.txt
OUTPUT_DIR=/data/rivera/lombard/out
mkdir -p ${OUTPUT_DIR}/fastqc_out

while read line_all
do
        SAMPLE_NAME=`echo ${line_all} | cut -d" " -f1`
        R1_FASTQ=`echo ${line_all} | cut -d" " -f2`
        R2_FASTQ=`echo ${line_all} | cut -d" " -f3`
	echo """
	module load fastqc
	fastqc -o ${OUTPUT_DIR}/fastqc_out ${R1_FASTQ} ${R2_FASTQ}
	""" > ../bsubFiles/fastqc_${SAMPLE_NAME}.bsub	
done < ${FILE_SAMPLE_LIST}
