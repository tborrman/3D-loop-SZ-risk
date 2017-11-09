#!/usr/bin/env python
import subprocess

for i in range(1, 18):
	subprocess.call('/data/zusers/borrmant/borrmant.ghpcc.project/Research/Schahram/Schahram-project/overlap_loops_eQTL.py ' +
	'-i clean_master_requested_loops_' + '{:02d}'.format(i) + ' -o overlap_loops_eQTL_' + str(i) + ' &', shell=True)
