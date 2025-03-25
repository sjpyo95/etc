import sys, os
import vcf

gatk_filter = ['orientation_bias', 'orientation', 't_lod', 'read_position', 'contamination', 'str_contraction', 'germline', 'panel_of_normals', 'position', 'strand_bias']

vartype_filter = ['synonymous_variant', 'intron_variant', 'intergenic_region', 'conservative_inframe_deletion', 'conservative_inframe_insertion', 'disruptive_inframe_insertion', 'disruptive_inframe_deletion']
def knownCHIP(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip()
		rows = [':'.join(x.strip().split('\t')) for x in f]
	return rows

def filter_somatic(vcf_file, out_file, knownchip):
	filts = []
	reader = vcf.Reader(open(vcf_file, 'r'))
	writer = vcf.Writer(open(out_file, 'w'), reader)
	knownchipN = 0; novelchipN = 0
	filtcount = 0
	for rec in reader:
		filt = rec.FILTER
#		if not filt in filts: filts += filt
##Filter out germline, PON, slippage flagged variants
		common_elmnt = list(set(filt) & set(gatk_filter))
		if len(common_elmnt) > 0 : filtcount += 1; continue
#		if 'germline' in filt: continue
#		elif 'panel_of_normals' in filt: continue
#		elif 'orientation_bias' in filt: continue
#		elif 'strand_artifact' in filt: continue
#		elif 'contamination' in filt: continue
		
#		elif 'slippage' in filt: continue
#		elif 'weak_evidence' in filt: continue

##Filter out variants with a allele frequency greater than 0.01 in gnomAD and NARD2 database
##Check if the variant is known CHIP driver variant
		chr = rec.CHROM; pos = str(rec.POS)
		chrpos = chr+':'+pos
		if chrpos in knownchip: 
			knownchipN += 1
			rec.INFO['CHIP'] = 'O'
			writer.write_record(rec)

##Find novel CHIP variants
		else:
			gnomad_all = rec.INFO['gnomAD_exome_ALL']
			gnomad_eas = rec.INFO['gnomAD_exome_EAS']
 			nard2 = rec.INFO['NARD2']
			exac_all = rec.INFO['ExAC_ALL']
#			af_max = 
			if exac_all is not None:
				if exac_all > 0.01 : continue
			if gnomad_all is not None:
				if gnomad_all > 0.01: continue
			elif gnomad_eas is not None:
				if gnomad_eas > 0.01: continue
#			elif nard2 is not None:
#				if nard2 > 0.01: continue

##Filter out the variants with a read depth (count) less than 10
			reads = []
			for sample in rec.samples:
				refn = sample['AD'][0]
				altn = sample['AD'][1]
				dp = sample['DP']
				if refn == 0 : vaf = 1
				else:
					vaf = float(altn)/refn
				if dp < 100 : continue
				##For CHIP variant
				if vaf < 0.01: continue
#				elif vaf < 0.01 : continue
#				elif vaf > 0.5: continue

#				ads = sample['AD']
#				if ads != None:
#					reads.append(ads[1])
#			read_avg = float(sum(reads))/len(reads)
#			if read_avg < 10: continue
			snpeff_anno = rec.INFO['ANN'][0]
			vartype = snpeff_anno.split('|')[1]
			if not vartype in filts: filts.append(vartype)
#			common_vartype = list(set(vartype) & set(vartype_filter))
			if vartype in vartype_filter: continue
#			if 'synonymous_variant' in vartype: continue
#			elif 'intron_variant' in vartype: continue
#			elif 'intergenic_region' in vartype: continue
#			elif '5_prime_UTR' in vartype: continue
#			elif 'conservative_inframe_deletion' in vartype: continue
			
##Filter out that does not have cosmic feature
#			cosmic = rec.INFO['cosmic70']
#			if None in cosmic: continue
#			print cosmic

##Filter out synonymous variants
			func = rec.INFO['ExonicFunc.refGene']
			if None in func: continue
#			elif 'weak_evidence' in filt: continue
#			print func
#			print filt
			
			novelchipN += 1
			rec.INFO['CHIP'] = 'X'
			writer.write_record(rec)

	writer.close()
#	print  knownchipN,'\t',novelchipN
#	return knownchipN, novelchipN
	return filts

inputdir = sys.argv[1]+'/'
outputdir = sys.argv[2]+'/'
knownfile = sys.argv[3]
knownchip = knownCHIP(knownfile)

if not os.path.exists(outputdir): os.makedirs(outputdir)

#vcf_files = [inputdir + x for x in filter(lambda x: '.vcf' in x and '.idx' not in x, os.listdir(inputdir))]
samples = filter(lambda x: '.vcf' not in x and 'header' not in x and not 'logs' in x, os.listdir(inputdir))
#for vcf_file in vcf_files:
#print 'Sample\tKnown\tNovel'
totalknown = 0; totalnovel = 0
#print 'Sample\tknownCHIP\tnovelCHIP'
filts_all = []
for sample in samples:
#	print sample+'\t',
	vcf_file = inputdir + sample + '/' + filter(lambda x: sample in x and '.vcf' in x and '.tbi' not in x, os.listdir(inputdir + sample + '/'))[0]
#	sample = vcf_file.split(inputdir)[-1].split('.m2')[0]
	soutputdir = outputdir + sample +'/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	out_file = soutputdir + sample + '.somatic.vcf'
	
#	knownN, novelN = filter_somatic(vcf_file, out_file, knownchip)
	filts = filter_somatic(vcf_file, out_file, knownchip)
	filts_all += filts
#	totalknown += knownN; totalnovel += novelN

#print list(set(filts_all))
#print 'Total\t',totalknown,'\t',totalnovel
