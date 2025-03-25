import sys, os
import re

def parse_panel(infile):
	group = {}
	with open(infile) as f:
		lines = f.readlines()
		for i in range(1, len(lines)):
			line = lines[i].strip()
			col = line.split('\t')
			samid = col[1]; sex = col[2]; age = col[3]
			idnum = int(id.split('SCA')[-1])
			samid = 'SCA-'+str(idnum)
			chcount = int(cols[4])
			if not group.has_key('CH'):
				group['CH'] = []
			if not group.has_key('Non-CH'):
				group['Non-CH'] = []
			if chcount > 0:
				group['CH'].append(samid)

			else:
				group['Non-CH'].append(samid)

	return group

def parse_array(infile):
	with open(infile) as f:
		lines = f.readlines()
		samples = lines[0].strip().split('\t')[1:]
		
		for i in range(1, len(lines)):
			line = lines[i].strip()
			col = line.split('\t')
			pos = col[0]; calls = col[1:]
			
