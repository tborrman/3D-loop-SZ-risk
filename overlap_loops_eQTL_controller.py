#!/usr/bin/env python
import subprocess

for i in range(18):
	subprocess.call('bsub -q long -W 12:00 -R "rusage[mem=5000]" -o bsub'+str(i)+'.out -e bsub'+str(i)+'.err ' +
		'/home/tb37w/project/Research/Schahram/Schahram-project/overlap_loops_eQTL.py ' +
		'-i clean_master_requested_loops_' + '{:02d}'.format(i) + 
		' -e cis_eQTL_PGC_SCZ2_filter_FDR_1e-20.tsv -o overlap_loops_eQTL_' + str(i), shell=True)
