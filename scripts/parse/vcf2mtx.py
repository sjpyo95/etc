import vcf
import sys
import csv

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


vcffile = sys.argv[1]
bedfile = sys.argv[2]
outfilename = sys.argv[3]

bed = parseBed(bedfile)

vcf_reader = vcf.Reader(open(vcffile, 'r'))
header = ['VAR_ID', 'GENE', 'VAR_PROP', 'AF']
samples = vcf_reader.samples
header.extend(samples)

data = []
for record in vcf_reader:
	if not record.CHROM in bed.keys(): gene = 'NA'
	else:
		region = filter(lambda x: record.POS >= x[0] and record.POS <= x[1], bed[record.CHROM])
		if len(region) == 0: gene = 'NA'
		else: gene = region[0][-1]
	if record.INFO.has_key('AF'): af = record.INFO['AF'][0]
	else: af = '0'
#	print record.CHROM
#	print record.POS
#	print record.REF
#	print record.ALT
	varid = record.CHROM+':'+str(record.POS)+':'+record.REF+'>'+str(record.ALT)
	row = [varid, gene, 0, af]
	for sample in record.samples:
		gt = sample['GT']
		ad = sample['AD']
		dp = sample['DP']
		af = float(ad[1])/float(dp)
#		row.append(gt)
		if gt != '0/0':
			row[2] += 1
			print af
			row.append(str(af*100))
#			row.append(sample.gt_bases[0]+sample.gt_bases[-1])
#			if gt == '0/1':
#				row.append(str(af))
#				row.append(sample.gt_bases[0] + sample.gt_bases[-1])
#			elif gt == '1/0':
#				row.append('1')
#				row.append(sample.gt_bases[0]+sample.gt_bases[-1]
#			elif gt == '1/1':
#				row.append('2')
#				row.append(sample.gt_bases[0]+sample.gt_bases[-1])
			
#			else:
#				row.append('.')
		else:
			row.append(str(af*100))
#			row.append('0')
	data.append(row)

#for row in data:
#	row[3] = str(row[3]) + '/' + str(len(samples))
outfile = open(outfilename, 'w')
outfile.write('\t'.join(header))
for row in data:
	outfile.write('\n'+ '\t'.join(map(str,row)))

