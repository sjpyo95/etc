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


header = ['VAR_ID', 'RS_ID', 'FILT', 'GENE_NM', 'CHIP', 'AA_CHANGE', 'Func', 'Pred', 'gnomAD_ALL', 'gnomAD_EAS', 'NARD2', 'COSMIC', 'CLNSIG', 'CLNDN']

inputdir = sys.argv[1]
samples = filter(lambda x: '.txt' not in x and '.vcf' not in x and 'log' not in x, os.listdir(inputdir))
samples = sorted(samples, key=lambda x: int(x.split('-')[-1]))
knownfile = sys.argv[2]
#knownfile2 = sys.argv[3]
outfile = sys.argv[3]

knownPos = knownCHIP(knownfile)
#knownAAchange = knownCHIP_AA(knownfile2)

header+=samples
newdic = {}
features = {}

for i in range(len(samples)):
	sample = samples[i]
	vcf_file = inputdir + sample + '/' + filter(lambda x: '.vcf' in x, os.listdir(inputdir+sample))[0]
	print sample
	reader = vcf.Reader(open(vcf_file, 'r'))
	for record in reader:
		chrom = record.CHROM; pos = record.POS
		ref = record.REF; alt = ','.join(list(map(str,record.ALT))); filt = record.FILTER; rsid = record.INFO['avsnp150'][0]
		varid = str(chrom)+':'+str(pos)+':'+ref+'>'+alt
		if rsid == None: rsid = '.'
		gene = record.INFO['Gene.refGene']
		gnomad_all = none2dot(record.INFO['gnomAD_exome_ALL'])
		gnomad_eas = none2dot(record.INFO['gnomAD_exome_EAS'])
		nard = none2dot(record.INFO['NARD2'][0])
		sift = none2dot(record.INFO['SIFT_pred'][0])
		fathmm = none2dot(record.INFO['FATHMM_pred'][0])
		provean = none2dot(record.INFO['PROVEAN_pred'][0])
		meta = none2dot(record.INFO['MetaSVM_pred'][0])
		mcap = none2dot(record.INFO['M-CAP_pred'][0])
		ai = none2dot(record.INFO['PrimateAI_pred'][0])
		deogen = none2dot(record.INFO['DEOGEN2_pred'][0])
		bayes = none2dot(record.INFO['BayesDel_addAF_pred'][0])
#		hdiv = none2dot(record.INFO['Polyphen2_HDIV_pred'][0])
#		hvar = none2dot(record.INFO['Polyphen2_HVAR_pred'][0])
#		lrt = none2dot(record.INFO['LRT_pred'][0])
#		pred = sift+';'+hdiv+';'+hvar+';'+lrt
		pred = list(set([sift, fathmm, provean, meta, mcap, ai, deogen, bayes]))
		if 'D' in pred: pred = 'D'
		elif '.' in pred: 
			if len(pred) == 1: pred = '.'
			else: pred.remove('.'); pred = ','.join(pred)
		else: pred = ','.join(pred)
		func = record.INFO['ExonicFunc.refGene']
#		print record.INFO['AAChange.refGene']

		if None in record.INFO['AAChange.refGene'] : aachanges = '.'
		else:
			aachanges = list(set(['p.'+x.split('p.')[-1] for x in record.INFO['AAChange.refGene']]))
#		print aachanges
		
		sample = record.samples[0].sample.split("-")[0]+'-'+record.samples[0].sample.split("_")[1]
		call = record.genotype(record.samples[0].sample)
		gt = call.data.GT
		ads = call.data.AD
		vaf = float(ads[1])/call.data.DP
#		print (chrom,pos,ref,alt)

		cosmic = record.INFO['cosmic70']
#		print cosmic
		if None in cosmic: cosmic = '.'
		else:
			cosmic = [x.split('(')[-1].replace(')','') for x in cosmic if '(' in x]
			cosmic = ','.join(cosmic)
#		print cosmic

		if not newdic.has_key((chrom,pos,ref,alt)):
			newdic[(chrom,pos,ref,alt)] = {x:'.' for x in samples}
		
		if len(filt) == 0: filt = ['PASS']
#		newdic[(chrom,pos,ref,alt)][sample] = (ads[0])+','+str(ads[1])+';VAF='+str(round(vaf,5))
		newdic[(chrom,pos,ref,alt)][sample] =str(round(vaf,5))
		
		if str(chrom)+':'+str(pos) in knownPos: 
			known = 'O'
#			if str(chrom)+':'+str(pos) in knownAAchange:
#				if bool(set(knownAAchange[str(chrom)+':'+str(pos)]) & set(aachanges)):
#					known = 'O:'+ ','.join(list(set(knownAAchange[str(chrom)+':'+str(pos)]) & set(aachanges)))
#				else: known = 'O'
#			else: known = 'O'	
		else: known = 'X'
#		print known

		if None in func:
			n = func.index(None)
			func[n] = '.'
		
		clnsig = record.INFO['CLNSG']
		clndn = record.INFO['CLNDN']
		if not features.has_key((chrom,pos,ref,alt)):
			features[(chrom,pos,ref,alt)] = [varid, rsid, filt, ','.join(gene), known, ','.join(aachanges),','.join(filter(None,func)), pred, gnomad_all, gnomad_eas, nard, cosmic, clnsig, clndn]
#		print filt
#		if 'PASS' in filt: features[(chrom,pos,ref,alt)][6] += 1
#		elif 'germline' in filt: features[(chrom,pos,ref,alt)][7] += 1
#		elif 'normal_of_panel' in filt: features[(chrom,pos,ref,alt)][8] += 1
#		else: features[(chrom,pos,ref,alt)][9] += 1

#		if gt != '0/0': features[(chrom,pos,alt)][8] += 
#		print features[(chrom,pos,ref,alt)]
#		print '\n'
chrposes = newdic.keys()
srt_pos = sorted(chrposes, key=lambda x: (int(x[0]) if x[0].isdigit() else float('inf'),int(x[1])))

lines = []
for chrpos in srt_pos:
	line = ''
	sampledic = newdic[chrpos]

	feature = list(map(str,features[chrpos]))
	line += '\t'.join(feature)
	for sample in samples:
		val = sampledic[sample]
		line += '\t' + val
	lines.append(line)

with open(outfile, 'w') as f:
	f.write('\t'.join(header)+'\n'+'\n'.join(lines))

#			print chrom, pos, newdic[chrom,pos]
