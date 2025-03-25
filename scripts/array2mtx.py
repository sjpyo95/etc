import sys, os
import csv
#sys.path.insert(0, '/home/sumyme95/scripts/modules/')
#import config as cf
import time
import gzip
import commands

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

def scoreCalls(call):
	if len(call) > 2: print call
#	elif call[0] == call[1]:
		
def parseBed(filename):
	beddic = {}
	with open(filename, 'rb') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			gene = row[3].split('.')[0]
			if not beddic.has_key(row[0]):
				beddic[row[0]] = []
			beddic[row[0]].append([int(row[1]), int(row[2]), gene])
	return beddic

def parseMacrogen(filename, genome, bed):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	regiondic = dict()
	refdic = dict()
	for i in range(len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		region_name = tmp[1]
		chr = tmp[9]; pos = int(tmp[10])
		if chr not in bed.keys(): continue
		targetRs = bed[chr]
		target = filter(lambda x: pos >= x[0] and pos <= x[1], targetRs)
		if len(target) == 0 : gene = 'NA'
		else: gene = target[0][-1]
		position = chr + ':' + str(pos)
		strand = tmp[2]
		regiondic[region_name] = position
		ref = genome[chr][int(pos)-1]
		refdic[region_name] = [ref, gene]
	return regiondic, refdic


def parseArray(filename):
	infile = open(filename, 'r')
	lines = infile.readlines()
	cols = lines[0].split('\t')
	samples = cols[1:]
	array = dict()
	newSams = []
	varcount = dict()
	for i in range(1,len(lines)):
		line = lines[i].strip()
		tmp = line.split('\t')
		pos = tmp[0]; calls = tmp[1:]
		varcount[pos] = 0
#		if not pos in refdic.keys(): continue
		
#		ref = refdic[pos]
		if not array.has_key(pos):
			array[pos] = dict()
		varcount[pos] = []
		for t in range(len(samples)):
			sample = samples[t].replace('-Chip', '').strip()
			newSams.append(sample)
			call = calls[t]
			array[pos][sample] = call
#			v = list(set(list(call)))
#			if call != ref+ref:
			varcount[pos].append(call)

#		varcount[pos] = list(set(list(varcount[pos])))
	return array, list(set(newSams)), varcount

def intersect(some_dict, another_dict):
	intersect = []
	for item in some_dict.keys(  ):
		if another_dict.has_key(item):
			intersect.append(item)
	return intersect

def sortChrs(chromosome_positions):
    return sorted(chromosome_positions, key=lambda x: (x[0].split(':')[0], int(x[0].split(':')[1])))

def writeMatrics(dic, outname, samples, positions, regiondic, refdic, varcount):
	outfile = open(outname, 'w')
	outfile.write('POS\tREF\tALT\tGENE\tVAR_PROP\t')
	outfile.write('\t'.join(samples)+'\n')
	lines = []
	for positionid in positions:
		line = ''
		samcalls = dic[positionid]
		position = regiondic[positionid]
		
		ref = refdic[positionid][0]

		gene = refdic[positionid][1]
		varlist = varcount[positionid]

#		print ref, len(filter(lambda x: x != ref+ref, varlist))

		line += position+'\t'+ref+'\t'
#		outfile.write(position+'\t'+ref+'\t')
		
		all_vars = list(set(list(''.join(samcalls.values()))))
		if ref in all_vars:
			all_vars.remove(ref)
		if len(all_vars) == 0:
#			outfile.write(ref)
			line += ref
		else:
#			outfile.write(','.join(all_vars))
			line += ','.join(all_vars)
		line += '\t'+gene
		count = len(filter(lambda x: x != ref+ref, varlist))
		line += '\t'+str(count)
#		outfile.write('\t' + str(count))
		for sample in samples:
			call = samcalls[sample]
#			outfile.write('\t'+str(call))
			
			line += '\t'+str(call)
#		line+= '\n'
		lines.append(line)
#		outfile.write('\n')
	outfile.write('\n'.join(lines))
	outfile.close()

if __name__ == '__main__' :

	arrayfile = sys.argv[1]
	reportfile = sys.argv[2]
	fastafile = sys.argv[3]
	bedfile = sys.argv[4]
	outfile = sys.argv[5]

#	sampleinfofile = sys.argv[3]
#	target = sys.argv[4]
	
	with open(fastafile) as fp:
		genome = dict(read_fasta(fp))
	
	bed = parseBed(bedfile)

	regioninfo, refdic = parseMacrogen(reportfile, genome, bed)

	print '\nGenome parsing fininshed.\n'
	array, samples, varcount = parseArray(arrayfile)
#	array, samples = parseArray(arrayfile)
	print 'Genome and Array file Parsing Finished!\n'
	samples.sort(key = lambda x: int(x.split('_')[-1]))
	
	int_regions = intersect(array, regioninfo)

#	filt_refdic = {k:v for k, v in refdic.items() if k in int_regions}
	print 'Writing Matrix...\n'
	writeMatrics(array, outfile, samples, int_regions, regioninfo, refdic, varcount)
