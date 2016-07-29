for f in  ../bsubFiles/hisat2_alignment_hg19_SRR*bsub
do
        jobname=`basename $f | sed 's/\.bsub//g'`
        cat $HOME/bsub_options $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
