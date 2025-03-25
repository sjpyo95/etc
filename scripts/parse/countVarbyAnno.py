import sys, os
import csv

mtxfile = sys.argv[1]
annfile = sys.argv[2]

with open(mtxfile, 'rb') as f:
	header = f.readline().strip().split('\t')
	mtx = [line.strip().split('\t') for line in f]


with open(annfile, 'rb') as f:
	lines = f.readlines()
	ann_dict = {}
	for i in range(len(lines)):
		line = lines[i].strip().split('\t')
		if not ann_dict.has_key(line[1]):
			ann_dict[line[1]] = []
		ann_dict[line[1]].append(line[0])

print mtx
samples = header[5:]
with open(outfile, 'w') as f:
	outfile.write('age\tcount\n')
	for age in ann_dict.keys():
		insams = ann_dict[age]
		for sam in insams:
			samindex = header.index(sam)
			
			
			
