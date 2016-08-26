for f in  ../bsubFiles/map_invalid_barcode1.bsub
do
        jobname=`basename $f | sed 's/\.bsub//g'`
        cat $HOME/bsub_options $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
