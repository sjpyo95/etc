import vcf
import sys, os

vcf_file = sys.argv[1]



header = ['CHROM', 'POS','rsID', 'REF', 'ALT', 'GENE', 'PASS','GERM','PON','ETC','CHIP', 'AA_CHANGE', 'Func', 'Pred', 'gnomAD_ALL', 'gnomAD_EAS', 'NARD2', 'COSMIC']
reader = vcf.Reader(open(vcf_file, 'r'))
for record in reader:
	chrom = record.CHROM; pos = record.POS
	ref = record.REF; alt = ','.join(list(map(str,record.ALT))); filt = record.FILTER; rsid = record.INFO['avsnp150'][0]
	if rsid == None: rsid = '.'
	gene = record.INFO['Gene.refGene']
	gnomad_all = none2dot(record.INFO['gnomAD_exome_ALL'])
	gnomad_eas = none2dot(record.INFO['gnomAD_exome_EAS'])
	nard = none2dot(record.INFO['ALL'][0])
	sift = none2dot(record.INFO['SIFT_pred'][0])
	fathmm = none2dot(record.INFO['FATHMM_pred'][0])
	provean = none2dot(record.INFO['PROVEAN_pred'][0])
	meta = none2dot(record.INFO['MetaSVM_pred'][0])
	mcap = none2dot(record.INFO['M-CAP_pred'][0])
	ai = none2dot(record.INFO['PrimateAI_pred'][0])
	deogen = none2dot(record.INFO['DEOGEN2_pred'][0])
	bayes = none2dot(record.INFO['BayesDel_addAF_pred'][0])
#	hdiv = none2dot(record.INFO['Polyphen2_HDIV_pred'][0])
#	hvar = none2dot(record.INFO['Polyphen2_HVAR_pred'][0])
#	lrt = none2dot(record.INFO['LRT_pred'][0])
#	pred = sift+';'+hdiv+';'+hvar+';'+lrt
	pred = list(set([sift, fathmm, provean, meta, mcap, ai, deogen, bayes]))
	if 'D' in pred: pred = 'D'
	elif '.' in pred: 
		if len(pred) == 1: pred = '.'
		else: pred.remove('.'); pred = ','.join(pred)
	else: pred = ','.join(pred)
	func = record.INFO['ExonicFunc.refGene']
#	print record.INFO['AAChange.refGene']

	if None in record.INFO['AAChange.refGene'] : aachanges = '.'
	else:
		aachanges = list(set(['p.'+x.split('p.')[-1] for x in record.INFO['AAChange.refGene']]))
#	print aachanges
	samples = record.samples
	for sample in samples:
		print sample

