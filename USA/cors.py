import glob

while len(glob.glob('/data/users/kgruber/results/USA/correlations/USA*.png'))<6:
	exec(open('/data/users/kgruber/Skripts/other/plot_USA_correlations.py').read())