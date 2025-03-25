import sys
import csv
def getAge(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		samAge = {}
		for line in f:
			samAge[line.strip().split('\t')[0]] = [line.strip().split('\t')[1],line.strip().split('\t')[2]]
	return samAge


parlist = sys.argv[1:]
mtxfile = parlist[0]; agefile = parlist[1]; outfile = parlist[2]

ages = getAge(agefile)

with open(mtxfile, 'r') as f:
	header  = f.readline().strip().split('\t')
	matrix = [line.strip().split('\t') for line in f]

samples = header[6:]
data = {}
data['CH_ALL'] = []
for row in matrix:
	gene = row[3]
	var = row[6:]
	if not data.has_key(gene):
		data[gene] = []
	for i in range(len(var)):
		v = var[i]; sample = samples[i]
		age = ages[sample][0]; sex = ages[sample][1]
		if v != '.':
			data[gene].append(age)
			data['CH_ALL'].append(age)
with open(outfile, 'w') as f:
	f.write('GENE\tAGE')

	for gene in data.keys():
		geneages = data[gene]
		for age in geneages:
			f.write('\n'+gene+'\t'+age+'\t')

