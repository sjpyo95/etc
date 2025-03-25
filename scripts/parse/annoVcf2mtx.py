import sys, os
import vcf

def none2dot(feature):
	if feature == None:
		feature = '.'
	return feature

def knownCHIP(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip()
		rows = [':'.join(x.strip().split('\t')) for x in f]
	return rows

header = ['VAR_ID', 'RS_ID','CHIP', 'sample_cnt', 'FILT', 'GENE_NM', 'VAF_MN', 'GERM_MX_AF', 'HGVS_CDNA', 'HGVS_PROTEIN', 'VAR_TYPE', 'Func', 'CLNSIG', 'CLNDN', 'COSMIC']
vcffile = sys.argv[1]
knownfile = sys.argv[2]
outfile = open(sys.argv[3], 'w')

reader = vcf.Reader(open(vcffile,'r'))
samples = reader.samples
#samples = sorted(samples, key=lambda x: int(x.split('-')[1]))
header.extend(samples)
knownPos = knownCHIP(knownfile)
outfile.write('\t'.join(header))
lines = ''
for record in reader:
	chr = record.CHROM; pos = str(record.POS)

	ref = record.REF; alt=','.join(list(map(str,record.ALT))); filt = record.FILTER
	if len(filt) == 0: filt.append('PASS')
#	if 'germline' in filt: continue
#	elif 'panel_of_normals' in filt: continue
#	elif 'slippage' in filt: continue

	rsid = record.INFO['avsnp150'][0]
	varid = str(chr)+':'+str(pos)+':'+ref+'>'+alt
	if rsid == None: rsid = '.'
	gene = ';'.join(record.INFO['Gene.refGene'])
	
	gnomad_all = none2dot(record.INFO['gnomAD_exome_ALL'][0])
	gnomad_eas = none2dot(record.INFO['gnomAD_exome_EAS'][0])
	nard = none2dot(record.INFO['NARD2'][0])
	exac = none2dot(record.INFO['ExAC_ALL'][0])
	germs = [s for s in [gnomad_all, gnomad_eas, nard, exac] if s != '.']
	if len(germs) == 0: germ_af = '.'
	else:
		germ_af = str(max(map(float, germs)))
	
	
#	exac = none2dot(record.INFO['ExAC_ALL'])
#	germ_af = str(gnomad_all)+';'+str(gnomad_eas)+';'+str(nard)

#	sift = none2dot(record.INFO['SIFT_pred'][0])
#	fathmm = none2dot(record.INFO['FATHMM_pred'][0])
#	provean = none2dot(record.INFO['PROVEAN_pred'][0])
#	meta = none2dot(record.INFO['MetaSVM_pred'][0])
#	mcap = none2dot(record.INFO['M-CAP_pred'][0])
#	ai = none2dot(record.INFO['PrimateAI_pred'][0])
#	deogen = none2dot(record.INFO['DEOGEN2_pred'][0])
#	bayes = none2dot(record.INFO['BayesDel_addAF_pred'][0])
##	hdiv = none2dot(record.INFO['Polyphen2_HDIV_pred'][0])
##	hvar = none2dot(record.INFO['Polyphen2_HVAR_pred'][0])
##	lrt = none2dot(record.INFO['LRT_pred'][0])
#	pred = meta+';'+mcap+';'+ai+';'+deogen+';'+bayes
#	pred = list(set([sift, fathmm, provean, meta, mcap, ai, deogen, bayes]))
#	if 'D' in pred: pred = 'D'
#	elif '.' in pred: 
#		if len(pred) == 1: pred = '.'
#		else: pred.remove('.'); pred = ','.join(pred)
#	else: pred = ','.join(pred)

	clnsig = none2dot(record.INFO['CLNSIG'][0])
	clndn = none2dot(record.INFO['CLNDN'][0])
	
	snpeff_ann = record.INFO['ANN'][0].split('|')

	vartype = snpeff_ann[1]
	hgvs_cdna = snpeff_ann[9]; hgvs_pr = snpeff_ann[10]

	func = record.INFO['ExonicFunc.refGene']
	if None in func:
		n = func.index(None)
		func[n] = '.'
	func = ';'.join(func)

#	if None in record.INFO['AAChange.refGene'] : aachanges = '.'
#	else:
#		aachanges = ';'.join(list(set(['p.'+x.split('p.')[-1] for x in record.INFO['AAChange.refGene']])))
	cosmic = record.INFO['cosmic70']
	if None in cosmic: cosmic = '.'
	else:
		cosmic = [x.split('(')[-1].replace(')','') for x in cosmic if '(' in x]
		cosmic = ','.join(cosmic)

	if str(chr)+':'+str(pos) in knownPos: 
		known = 'O'
	else: known = 'X'
	
	sample_cnt = 0
	vafs = []
	vafs4mean = []
	for sample in record.samples:
		ads = sample['AD']
		if ads == None: 
			vaf = '.'
		else:
			refad = float(ads[0]); varad = float(ads[1])
			vaf = float(varad)/(refad+varad)
			sample_cnt+=1
#			if vaf >= 0.02: sample_cnt+=1
			vaf = round(vaf,4)
			vafs4mean.append(vaf)

		vafs.append(str(vaf))
		
	vafmean = round(sum(vafs4mean)/len(vafs4mean), 4)
#	print varid, rsid, str(sample_cnt), ';'.join(filt), gene, str(vafmean), known, hgvs_cdna, hgvs_pr, vartype, func, germ_af, cosmic, clnsig, clndn, '\t'.join(vafs)
	line = '\t'.join([varid, rsid, known, str(sample_cnt), ','.join(filt), gene, str(vafmean), germ_af, hgvs_cdna, hgvs_pr, vartype, func, clnsig, clndn, cosmic, '\t'.join(vafs)])
	lines+='\n'+line

outfile.write(lines)

outfile.close()

	





