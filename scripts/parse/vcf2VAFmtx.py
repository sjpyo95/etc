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

#def sorting(l):
#	return sorted(l, key = lambda x: x[-1])
#filename = sys.argv[1]
#
#reader = vcf.Reader(open(filename, 'r'))
#
#het = []; hom = []; ref = []
#for record in reader:
#	for sample in record.samples:
#		ad = sample['AD']
#		dp = sample['DP']
#		gt = sample['GT']
##		vaf = float(ad[gt)/dp
#		if gt != '0/0' and gt != './.':
#			if len(ad) > 2: continue
#				
#				#hom.append([gt,vaf])
#			elif gt == '0/1' or gt == '1/0':
#				vaf = float(ad[0])/dp
#				het.append([sample,gt,ad, dp,vaf, record])
#			else:
#				vaf = float(ad[0])/dp
#				hom.append([sample,gt,ad,dp,vaf, record])
#		
#		else: pass
#			ref.append([gt,vaf])
#
#sort_het =  sorting(het)
#for i in sort_het:

filename = sys.argv[1]
bedfile = sys.argv[2]
outfile = sys.argv[3]
bed = parseBed(bedfile)
vcf_reader = vcf.Reader(open(filename, 'r'))
samples = vcf_reader.samples
#selsamples = ['SCA-IDT_109','SCA-IDT_238','SCA-IDT_249','SCA-IDT_286','SCA-IDT_295','SCA-IDT_247','SCA-IDT_284','SCA-IDT_308','SCA-IDT_260','SCA-IDT_241','SCA-IDT_182','SCA-IDT_248','SCA-IDT_298','SCA-IDT_181','SCA-IDT_307','SCA-IDT_237','SCA-IDT_306','SCA-IDT_254','SCA-IDT_253','SCA-IDT_304']
#samples_new = [x.split('_S')[0] for x in samples if x.split('_S')[0] in selsamples]
#samples = filter(lambda x: x.split('_S')[0] in selsamples, samples)
header = ['VAR_ID', 'GENE', 'VAR_PROP', 'AF']
header.extend(samples)
#header.extend(samples_new)
#print len(samples_new)
variants = []
het = []; hom = []
for record in vcf_reader:
	if not record.CHROM in bed.keys(): gene = 'NA'
	else:
		region = filter(lambda x: record.POS >= x[0] and record.POS <= x[1], bed[record.CHROM])
		if len(region) == 0: gene = 'NA'
		else: gene = region[0][-1]
	if record.INFO.has_key('AF') : af = record.INFO['AF'][0]
	else : af = '0'
	row = [record.CHROM + ':' + str(record.POS)+':'+record.REF+'>'+str(record.ALT[0]), gene, 0, af]
	for sample in samples:
		call = record.genotype(sample)
		if call['GT'] == './.':
			vaf = '.'
#		if call.gt_type == 0:
#			vaf = 0.0
#		elif call.gt_type == 2:
#			vaf = 1.0
		else:
			row[2] += 1
			af1 = float(call.data.AD[1])/call.data.DP
			af2 = float(call.data.AD[0])/call.data.DP
#			vaf = min([af1,af2])
			vaf = float(call.data.AD[1])/call.data.DP
			if call.gt_type == 1:
				het.append(vaf)
				pass
#				print call
			elif call.gt_type == 2:
				hom.append(vaf)
#				pass
#				print call
			
				
		row.append(vaf)
	variants.append(row)

#print max(het), min(het)
#print max(hom), min(hom)
with open(outfile, 'w') as f:
	f.write('\t'.join(header))
	for row in variants:
		f.write( '\n' + '\t'.join(map(str, row)))
f.close()
