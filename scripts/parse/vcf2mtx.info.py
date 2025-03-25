import vcf
import sys
import csv

def median(lst):
	n = len(lst)
	s = sorted(lst)
	if n % 2 == 0:
		return (s[n//2-1] + s[n//2])/2.0
	else:
		return s[n//2]

def getAnnovar(annovar):
	reader = vcf.Reader(open(annovar,'r'))
	data = {}
	for rec in reader:
		chrom = rec.CHROM; pos = rec.POS
		exacall = rec.INFO['ExAC_ALL']
		exaceas = rec.INFO['ExAC_EAS']
		nard = rec.INFO['ALL'][0]
		sift = rec.INFO['SIFT_pred'][0]
		hdiv = rec.INFO['Polyphen2_HDIV_pred'][0]
		hvar = rec.INFO['Polyphen2_HVAR_pred'][0]
		lrt = rec.INFO['LRT_pred'][0]
		data[chrom+':'+str(pos)] = [exacall, exaceas, nard, sift, hdiv, hvar, lrt]
	return data

def getBank(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		rows = [line.strip().split('\t') for line in f]
	data = {}
	for row in rows:
		chrom = row[0]; pos = row[1]
		gene = row[5]
		ref = row[2]; alt = row[3]
		aachange = row[-1]
		conseq = row[-2]
		chrom_pos = chrom + ':' + pos
		if not data.has_key(chrom_pos):
			data[chrom_pos] = []
		data[chrom_pos].append([ref,alt,aachange,conseq])
	return data

def getMaf(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		rows = [line.strip().split('\t') for line in f]
	data = {}
	for row in rows:
		chrom = row[4]; pos = row[5]
		conseq = row[8]
		data[chrom+':'+pos] = conseq
	return data

vcffile = sys.argv[1]
annovarfile = sys.argv[2]
maffile =sys.argv[3]
bankfile = sys.argv[4]
outfilename = sys.argv[5]

annovar = getAnnovar(annovarfile)
maf = getMaf(maffile)
bank = getBank(bankfile)

vcf_reader = vcf.Reader(open(vcffile, 'r'))
header = ['POS', 'REF', 'ALT', 'GENE', 'VAR_PROP', 'AF', 'VAF', 'CONSEQUENCE', 'KNOWN', 'ExonicFunc', 'Pred', 'ExAC_ALL', 'ExAC_EAS', 'NARD2']
samples = vcf_reader.samples
header.extend(samples)

data = []
for record in vcf_reader:
	chrom = record.CHROM; pos = record.POS
	anno = record.INFO['ANN'][0].split('|')
	gene = anno[3]
	func = anno[1]
	vafs = [float(record.genotype(x).data.AD[1])/record.genotype(x).data.DP for x in samples]
	vafs = filter(lambda x: x != 0, vafs)
	chr_pos = chrom+':'+str(pos)
	if chr_pos in annovar.keys():
		info = annovar[chr_pos]
		exacall = info[0]; exaceas = info[1]; nard = info[2]; pred = list(set(info[3:]))
	else: print annovar[chr_pos]
	if chr_pos in maf.keys():
		conseq = maf[chr_pos]
	else:
		conseq = '.'

#	if not record.CHROM in bed.keys(): gene = 'NA'
#	else:
#		region = filter(lambda x: record.POS >= x[0] and record.POS <= x[1], bed[record.CHROM])
#		if len(region) == 0: gene = 'NA'
#		else: gene = region[0][-1]
	
	if chr_pos in bank.keys():
		known = 'O'
		print record
		print bank[chr_pos]
	else:
		known = 'X'
	if record.INFO.has_key('AF'): af = record.INFO['AF'][0]
	else: af = '0'

	row = [record.CHROM + ':' + str(record.POS), record.REF, record.ALT[0], gene, 0, af, median(vafs), conseq, known, func, pred, exacall, exaceas, nard]

	for sample in record.samples:
		gt = sample['GT']
		if gt != '0/0':
			row[4] += 1
			row.append(sample.gt_bases[0]+sample.gt_bases[-1])
#			if gt == '0/1':
#				row.append(sample.gt_bases[0] + sample.gt_bases[-1])
#			elif gt == '1/0':
#				row.append(sample.gt_bases[0]+sample.gt_bases[-1]
#			elif gt == '1/1':
#				row.append(sample.gt_bases[0]+sample.gt_bases[-1])
#			
#			else:
#				row.append('.')
		else:
			row.append('.')
	
	data.append(row)

#for row in data:
#	row[3] = str(row[3]) + '/' + str(len(samples))
outfile = open(outfilename, 'w')
outfile.write('\t'.join(header))
for row in data:
	outfile.write('\n'+ '\t'.join(map(str,row)))

