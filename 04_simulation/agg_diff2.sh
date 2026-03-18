for thr in 20 30 40 50; do 
python agg_diff3.py ctcf_0pct wild 'CTCF(0%)-WT' $thr
python agg_diff3.py loading_0pct wild 'Cohesin(0%)-WT' $thr
python agg_diff3.py ctcf_0pct loading_0pct 'CTCF(0%)-Cohesin(0%)' $thr
done
