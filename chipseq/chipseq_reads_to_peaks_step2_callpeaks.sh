SAMPLE_NAME=$1
TREATMENT_FILE=$2
CONTROL_FILE=$3
OUTPUT_DIR=$4


module load python/2.7.3
module load macs2/2.1

macs2 callpeak -g hs -q 0.01 -t /PHShome/gz976/rivera/Data/08.11.16.Nextseq/MSC1002.KO.BAF47FLI1.FLI1.c.bam -c /PHShome/gz976/rivera/Data/08.11.16.Nextseq/MSC1002.KO.BAF47FLI1.INPUT.c.bam -n MSC1002.KO.BAF47FLI1.FLI1.c_VS_MSC1002.KO.BAF47FLI1.INPUT.c.q01 --outdir /PHShome/si992/baf_ews/results/final/macs2_out
macs2 callpeak -g hs -q 0.01 -t /PHShome/gz976/rivera/Data/08.11.16.Nextseq/MSC1002.KO.empty.FLI1.c.bam -c /PHShome/gz976/rivera/Data/08.11.16.Nextseq/MSC1002.KO.empty.INPUT.c.bam /PHShome/gz976/rivera/Data/08.11.16.Nextseq/MSC1002.KO.empty.INPUT.b.bam -n MSC1002.KO.empty.FLI1.c_VS_MSC1002.KO.empty.INPUT.b_plus_c.q01 --outdir /PHShome/si992/baf_ews/results/final/macs2_out
