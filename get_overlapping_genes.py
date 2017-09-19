#!/usr/bin/env python
import argparse
import sys
import os

# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser=argparse.ArgumentParser(description='Get genes overlapping loop loci for input cell type')
parser.add_argument('-g', help= 'Gene table (ex. gene_table)')
parser.add_argument('-l', help='Loop table (ex. clean_master_requested_loops)')
parser.add_argument('-c', help='Cell type (Astro, GM, Neu, or NPC )')
args = parser.parse_args()

def counter(i):
	if i%100 == 0:
		print 'On line: '+str(i)

def overlap(a, b):
        '''
        Return boolean for whether interval b
        overlaps interval a
        '''
        if (a[0] <= b[0] <= a[1]) or (b[0] <= a[0] <= b[1]):
                return True
        else:
                return False

def check_start_below_stop(start, stop):
	if start >= stop:
		print 'ERROR: start >= stop'
		sys.exit()

def is_loop_in_cell(line, cell):
	col_id = {'Astro': 21,
				'GM': 22,
				'Neu': 23,
				'NPC': 24}
	splitline = line.rstrip().split('\t')
	qval = float(splitline[col_id[cell]])
	if qval < -1:
		return True
	else:
		return False

def write_overlap_genes(cell, geneFile, loop_chrom, loop_interval, outfile, gene_set):
	'''
	Write genes overlapping loop calls
	to file
	'''
	with open(geneFile, 'r') as g:
		g.readline()
		for line in g:
			splitline = line.split('\t')
			ENSG = splitline[4]
			if ENSG not in gene_set:
				gene_chrom = splitline[0]
				gene_start = int(splitline[1])
				gene_stop = int(splitline[2])
				gene_interval = [gene_start, gene_stop]
				check_start_below_stop(gene_start, gene_stop)
				# Check if gene overlaps loop region
				if (gene_chrom == loop_chrom) and overlap(loop_interval, gene_interval):
					gene_set.add(ENSG)
					outfile.write(line)
	return gene_set
				

def main():
	# Initialize empty gene list
	genes = set()
	# OUTPUT file
	OUT = open(args.c + '_overlapping_genes.txt', 'w' )
	with open(args.l, 'r') as l:
		# Skip header
		l.readline()
		# Parse loop call table
		for i, loop_line in enumerate(l):
			counter(i)
			# Check if loop exists in cell type
			if is_loop_in_cell(loop_line, args.c):
				split_loop_line = loop_line.split('\t')
				loop_chrom = split_loop_line[0]
				loop_start = int(split_loop_line[1])
				loop_stop = int(split_loop_line[5])
				check_start_below_stop(loop_start, loop_stop)
				loop_interval = [loop_start, loop_stop]
				# Parse gene table for overlapping genes
				genes = write_overlap_genes(args.c, args.g, loop_chrom, loop_interval, OUT, genes)
				
	OUT.close()
if __name__ == '__main__':
	main()
