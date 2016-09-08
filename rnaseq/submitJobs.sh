for f in  ../bsubFiles/generatebw_*bsub
do
        jobname=`basename $f | sed 's/\.bsub//g'`
        cat $HOME/bsub_options $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
