#!/bin/bash

COMBINE_FASTQ="NO"
for i in "$@"
do
case $i in
    -g=*|--genome=*)
    GENOME="${i#*=}" #Zv9
    shift # past argument=value
    ;;
    -s=*|--sampleFile=*)
    FILE_SAMPLE_LIST="${i#*=}" #test_Zv9_config.txt
    shift # past argument=value
    ;;
    -o=*|--outputDir=*)
    OUTPUT_DIR="${i#*=}"  # /data/langenau/out
    shift # past argument=value
    ;;
    --combineFastq)
    COMBINE_FASTQ="YES"
    shift # past argument with no value
    ;;
    *)
            # unknown option
    ;;
esac
done

echo $GENOME $FILE_SAMPLE_LIST $OUTPUT_DIR $COMBINE_FASTQ
