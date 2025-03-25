import os, sys
import re
import timeit

arrayfile = sys.argv[1]
infofile = sys.argv[2]
fastafile = sys.argv[3]
outfile = sys.argv[4]

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name[1:].split(' ')[0], ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name[1:].split(' ')[0], ''.join(seq))

with open(fastafile) as fp:
	seqs = dict(read_fasta(fp))
print 'FASTA  Finished'

def getpos(posid,seqs):
	with open(infofile,'r') as f:
		lines = f.readlines()
		for i in range(len(lines)):
			line = lines[i].strip().split('\t')
			if posid == line[1]:
				chr = line[9]; pos = line[10]
				if chr not in seqs.keys(): continue
				ref = seqs[chr][int(pos)-1]
				pos = chr+':'+pos
				break
			else: pos = ''; ref = ''
	return pos,ref

def removeDup(str):
	s = set(str)
	s = ''.join(s)
	return s

def parseArrayResult(infile):
	array = {}
	with open(infile, 'r') as f:
		lines = f.readlines()
		names = lines[0].strip().split('\t')[1:]
		names = map(lambda x: re.sub(r'SCA-Chip_',r'SCA-', x), names)
		for i in range(1,len(lines)):
			print i,''
			line = lines[i].strip()
			col = line.split('\t')
			posid = col[0]; calls = col[1:]
			pos,ref = getpos(posid,seqs)
			alts = ','.join(set(''.join(calls)))

			if pos == '': continue

			if not array.has_key(pos):
				array[pos] = {}	
			for j in range(len(names)):
				name = names[j]
				call = calls[j]
				array[pos][name] = call
	timeit.timeit("parseArray", setup="from __main__ import parseArray, data", number=100)
	return array,names

array,names = parseArrayResult(arrayfile)

with open(outfile,'w') as f:
	f.write('ID\t'+'\t'.join(names)+'\n')
	for pos in array.keys():
		calls = array[pos]
		f.write(pos)
		for name in names:
			call = calls[name]
			f.write('\t'+call)
		f.write('\n')
