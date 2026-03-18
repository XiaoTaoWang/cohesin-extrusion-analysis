#python eg_diff.zoom.py ctcf_0pct wild 'CTCF(0%)-WT'
#python eg_diff.zoom.py loading_0pct wild 'Cohesin(0%)-WT'
#python eg_diff.zoom.py ctcf_0pct loading_0pct 'CTCF(0%)-Cohesin(0%)'
#python eg_diff.zoom.py ctcf_25pct wild 'CTCF(25%)-WT'
#python eg_diff.zoom.py ctcf_50pct wild 'CTCF(50%)-WT'
#python eg_diff.zoom.py ctcf_10pct wild 'CTCF(10%)-WT'

for treat in ctcf_10pct ctcf_25pct ctcf_50pct loading_10pct loading_25pct loading_50pct; do 
    python eg_diff.zoom2.py $treat wild 'CTCF(0%)-WT'
done
