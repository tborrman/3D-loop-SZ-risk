#!/usr/bin/env python
import argparse
import sys
parser = argparse.ArgumentParser(description='combine master_loops file with requested loop lists')
parser.add_argument('-i', help='input file (clean_master_requested_loops)', type = str, required=True)
args = parser.parse_args()


def get_resolution(line):
	start = int(line.split('\t')[1])
	end = int(line.split('\t')[2])
	r = end - start
	if r != 5000 and r != 10000:
		print 'ERROR: not 5kb or 10kb resolution'
		sys.exit()
	return r



def add_qvals(line, F):
	res = get_resolution(line)
	# Coordinates [chr1, x1, x2, chr2, y1, y2]
	coord = line.split('\t')[:6]
	values = []
	for cell in ['dopa_NurrPOS_NeuNPOS_9k', 'glia_Nurr1aneg_NeuNneg_9k', 'neurons_acc_9k']:
		found = False
		REQ = open('/home/tb37w/project/Research/Schahram/test/opt/juicer/work/sergio/downsampled/' + cell + '/hiccups/requested_list_' + str(res), 'r')
		for reqline in REQ:
			reqcoord = reqline.split('\t')[:6]
			if coord == reqcoord:
				found = True
				values = values + map(float,reqline.split('\t')[7:])
				break
		if not found:
			values = values + ["NA"]*13
		REQ.close()
	if len(values)!= 39:
		print 'ERROR: something went wrong with grabbing qvals'
		sys.exit()
	else:
		F.write('\t'.join(line.strip().split('\t')[:21] + map(str,values)) + '\n')

def main():
	OUT = open('sergio_requested_loops.txt', 'w')
	IN = open(args.i, 'r')
	header = IN.readline().strip().split('\t')[:21]
	RH = open('/home/tb37w/project/Research/Schahram/test/opt/juicer/work/sergio/downsampled/dopa_NurrPOS_NeuNPOS_9k/hiccups/requested_list_10000_oldlabel', 'r')
	Rheader = RH.readline()
	rheadlist = Rheader.strip().split('\t')[7:]
	add_header = []
	for cell in ['dopa_NurrPOS_NeuNPOS_9k', 'glia_Nurr1aneg_NeuNneg_9k', 'neurons_acc_9k']:
		for rh in rheadlist:
			add_header.append(cell + '_' + rh)
	OUT.write('\t'.join(header) + '\t' + '\t'.join(add_header) + '\n')
	for i, line in enumerate(IN):
		if i % 100 == 0:
			print 'On row: ' + str(i)
		add_qvals(line, OUT)

	IN.close()
	OUT.close()

if __name__ == '__main__':
	main()