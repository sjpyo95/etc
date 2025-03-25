import sys
import csv


def rsID(dbsnpfile):
	rsid_pos = {}
	with open(dbsnpfile, 'r') as f:
		header = f.readline().strip().split('\t')
		for line in f:
			tmp = line.strip().split('\t')
			rsid_pos[tmp[2]] = [tmp[0],tmp[1]]
	return rsid_pos

mtxfile = sys.argv[1]
rsid_pos = sys.argv[2]

with open(mtxfile, 'r') as f:
	header = f.readline().strip().split('\t')
	matrix = [line.strip().split('\t') for line in f]


samples = header[8:]

for row in matrix:
	

