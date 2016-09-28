for f in  ../bsubFiles/tcc_kallisto_*bsub
do
        jobname=`basename $f | sed 's/\.bsub//g'`
        cat $HOME/bsub_options_quick $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
