import sys
import vcf
import sys
import csv

snpeffvcf = sys.argv[1]
annovarvcf = sys.argv[2]
outfile = sys.argv[3]

with open('/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/1/somatic/targetseq_age.txt', 'r') as f:
	age_dict = dict(line.strip().split('\t') for line in f)

newData = []
reader = vcf.Reader(open(annovarvcf, 'r'))
synCount = dict()
dtCount = dict()
for rec in reader:
	chrom = rec.CHROM; pos = rec.POS
	exacall = rec.INFO['ExAC_ALL']
	exaceas = rec.INFO['ExAC_EAS']
	nard = rec.INFO['ALL'][0]
	af = max([float(i) for i in rec.INFO['AF']])
	syn = rec.INFO['ExonicFunc.refGene'][0]
#	print syn
	if not synCount.has_key(syn):
		synCount[syn] = 0
	dt = rec.INFO['SIFT_pred'][0]
	hdiv = rec.INFO['Polyphen2_HDIV_pred'][0]
	hvar = rec.INFO['Polyphen2_HVAR_pred'][0]
	lrt = rec.INFO['LRT_pred'][0]
	pred = [dt,hdiv,hvar,lrt]
#	polyphen = rec.INFO['polyphen2_HDIV_pred'][0]

	if not dtCount.has_key(dt):
		dtCount[dt] = 0
#	print rec, exacall, exaceas, nard, af
	if exacall is not None and exaceas is not None and nard is not None:
#		print rec, exacall, exaceas, nard
		if float(exacall) < 0.1 and float(exaceas) < 0.1 and float(nard) < 0.1:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
#			print rec, exacall, exaceas, nard
	elif exacall is not None and exaceas is not None:
		if float(exacall) < 0.1 and float(exaceas) < 0.1:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
#			print rec, exacall, exaceas, nard
	elif exacall is not None and nard is not None:
		if float(exacall) < 0.1 and float(nard) < 0.1:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
	elif exaceas is not None and nard is not None:
		if float(exaceas) < 0.1 and float(nard) < 0.1:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
	elif nard is not None:
		if float(nard) < 0.05:
#			print rec, exacall, exaceas, nard
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
	elif exaceas is not None:
		if float(exaceas) < 0.05:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
#			print rec, exacall, exaceas, nard
	elif exacall is not None:
		if float(exacall) < 0.05:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
#		else:
#			print rec, exacall, exaceas, nard
	else:
		if 'D' in pred or 'P' in pred:
			newData.append([chrom,pos])
			synCount[syn] += 1
			dtCount[dt] += 1
#			print rec, af
			pass
print '\n'
print synCount
print dtCount
print '\n'

snpReader = vcf.Reader(open(snpeffvcf, 'r'))
oldsamples = snpReader.samples
samples_sorted = sorted(snpReader.samples, key = lambda x: int(age_dict[x]))
snpReader.samples = samples_sorted
print snpReader.samples
writer = vcf.Writer(open(outfile, 'w'), snpReader)
snpReader.samples = oldsamples
for rec in snpReader:
	chrom_pos = [rec.CHROM, rec.POS]
	if chrom_pos not in newData: continue
	rec.samples = sorted(rec.samples, key = lambda x: int(age_dict[x.sample]))
	writer.write_record(rec)
writer.close()
			


