import vcf
import sys, os

def knownCHIP(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip()
		rows = [':'.join(x.strip().split('\t')) for x in f]
	return rows

def knownCHIP_AA(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		rows = {}
		for row in f:
			line = row.strip().split('\t')
			chrom = line[0]; pos = line[1]; aachange = line[-1]
			if not rows.has_key(chrom+':'+pos):
				rows[chrom+':'+pos] = []
			rows[chrom+':'+pos].append(aachange)
	return rows

def none2dot(feature):
	if feature == None:
		feature = '.'
	return feature

inputdir = sys.argv[1]

knownfile = sys.argv[2]
knownsites = knownCHIP(knownfile)

samples = filter(lambda x: '.vcf' not in x, os.listdir(inputdir))

print 'Sample\trsID\tCHROM\tPOS\tTYPE\tGENE\tCHIP\tGT\tVAR\tREF\tALT\tVAF\tgnomAD_ALL\tgnomAD_EAS\tFILT\tFUNC\tPRED\tCOSMIC'

for i in range(len(samples)):
	sample = samples[i]
	vcf_file = inputdir +'/'+ sample+ '/' + filter(lambda x: '.vcf' in x, os.listdir(inputdir+'/'+sample))[0]
	reader = vcf.Reader(open(vcf_file, 'r'))
	for record in reader:
		chrom = record.CHROM; pos = record.POS
		ref = record.REF; alt = record.ALT
		gene = record.INFO['Gene.refGene']
		mutation = [ref+'>'+str(x) for x in alt]
		mutation = ';'.join(mutation)
		if len(mutation) > 3: vtype = 'INDEL'
		else: vtype = 'SNV'

		rsid = record.INFO['avsnp150'][0]
		if rsid == None: rsid = '.'
		sample = record.samples[0].sample.split("-")[0]+'-'+record.samples[0].sample.split("_")[1]
		call = record.genotype(record.samples[0].sample)
		gt = call.data.GT
		ads = call.data.AD
		dp = call.data.DP
		vaf = round(float(ads[1])/dp, 5)
		filt = record.FILTER
		if len(filt) == 0: filt = ['PASS']
		if str(chrom)+':'+str(pos) in knownsites: known = 'O'
		else: known = 'X'
#gnomAD
		gnomad_all = none2dot(record.INFO['gnomAD_exome_ALL'])
		gnomad_eas = none2dot(record.INFO['gnomAD_exome_EAS'])
		func = record.INFO['ExonicFunc.refGene']
		if None in func:
			n = func.index(None)
			func[n] = '.'
		func = ','.join(filter(None,func))
#pred
		sift = none2dot(record.INFO['SIFT_pred'][0])
		fathmm = none2dot(record.INFO['FATHMM_pred'][0])
		provean = none2dot(record.INFO['PROVEAN_pred'][0])
		meta = none2dot(record.INFO['MetaSVM_pred'][0])
		mcap = none2dot(record.INFO['M-CAP_pred'][0])
		ai = none2dot(record.INFO['PrimateAI_pred'][0])
		deogen = none2dot(record.INFO['DEOGEN2_pred'][0])
		bayes = none2dot(record.INFO['BayesDel_addAF_pred'][0])
		pred = list(set([sift, fathmm, provean, meta, mcap, ai, deogen, bayes]))
		if 'D' in pred: pred = 'D'
		elif '.' in pred:
			if len(pred) == 1: pred = '.'
			else: pred.remove('.'); pred = ','.join(pred)
		else: pred = ','.join(pred)
#cosmic		
		cosmic = record.INFO['cosmic70']
		if None in cosmic: cosmic = '.'
		else:
			cosmic = [x.split('(')[-1].replace(')','') for x in cosmic if '(' in x]
			cosmic = ','.join(cosmic)

		print '\t'.join([sample, rsid, str(chrom), str(pos), vtype,','.join(gene), known, gt, mutation, str(ads[0]), str(ads[1]), str(vaf), str(gnomad_all), str(gnomad_eas), ','.join(filt), func, pred, cosmic])
