import sys, os
import re
import math

def parse_anno(infile):
	refdic = {}
	with open(infile) as f:
		lines = f.readlines()
		for i in range(len(lines)):
			line = lines[i].strip()
			if not line.startswith('#'):
				col = line.split('\t')
				chr = re.sub(r'chr', r'', col[0]); type = col[2]
				if not type == 'gene': continue
				start = int(col[3]); end = int(col[4])
				info = col[8]
				gene = info.split('gene_name "')[1].split('"')[0]
				if not refdic.has_key(chr):
					refdic[chr] = {}
				if not refdic[chr].has_key(gene):
					refdic[chr][gene] = []
				refdic[chr][gene].append([start,end])
	return refdic

def parse_samplefile(infile):
	group = {}
	with open(infile) as f:
		lines = f.readlines()
		for i in range(1, len(lines)):
			line = lines[i].strip()
			cols = line.split('\t')
			id = cols[1]; sex = cols[2]; age = cols[3]
			chcount = int(cols[4]); maxvaf = cols[5]
			id = re.sub(r'SCA', r'SCA-', id)
#			idnum = int(id.split('-')[-1])
#			id = 'SCA-'+str(idnum)
			if not group.has_key('CH'):
				group['CH'] = []
			if not group.has_key('Non-CH'):
				group['Non-CH'] = []
			if chcount > 0:
				group['CH'].append(id)
			else:
				group['Non-CH'].append(id)
	return group

def parse_vaf_mtx(mtxfile):
	mtx = {}
	with open(mtxfile) as f:
		lines = f.readlines()
		samples = lines[0].strip().split('\t')[4:]
		for i in range(1,len(lines)):
			line = lines[i].strip()
			cols = line.split('\t')
			var_id = cols[0]; rs_id = cols[1]
			chip = cols[2]
			vafs = cols[4:]
			if not mtx.has_key((var_id, chip)):
				mtx[(var_id,chip)] = {}
			for j in range(len(samples)):
				sample = samples[j]
				vaf = vafs[j]
				mtx[(var_id,chip)][sample] = vaf
	return mtx

def target_seq():
	vaf_mtx = open(outputdir+'/vaf.mtx.txt', 'w')
	fc_mtx = open(outputdir+'/fc.mtx.txt', 'w')
	vlines = []
	flines = []
	vaf_mtx.write('var_id\tCH\tNon-CH\t'+'\t'.join(ch_sam+non_ch_sam)+'\n')
	fc_mtx.write('var_id\tCH\tNon-CH\tCHavg\tNon-CHavg\tlog2FC\n')
	sam_ids_order = ch_sam+non_ch_sam
	for info in mtx.keys():
		samdic = mtx[info]
		var_id = info[0]; chip = info[1]
		vrow = [var_id, 0,0]
		for sam in sam_ids_order:
			vaf = samdic[sam]
			if vaf != '.':
				if sam in ch_sam:
					vrow[1] += 1
				else:
					vrow[2] += 1

			vrow.append(str(vaf))
		vlines.append(vrow)
		
		frow = [var_id, 0,0]
		ch_vaf_list = []
		for sam in ch_sam:
			vaf = samdic[sam]
			if vaf != '.':
				frow[1] += 1
				ch_vaf_list.append(float(vaf))

		nonch_vaf_list = []
		for sam in non_ch_sam:
			vaf = samdic[sam]
			if vaf != '.':
				frow[2] += 1
				nonch_vaf_list.append(float(vaf))
	
		if len(ch_vaf_list) > 0:
			ch_vaf_avg = sum(ch_vaf_list)/len(ch_vaf_list)

		else: ch_vaf_avg = 0; fc= 0

		if len(nonch_vaf_list) > 0:
			non_vaf_avg = sum(nonch_vaf_list)/len(nonch_vaf_list)
		else: non_vaf_avg = 0; fc = 'inf'
		if ch_vaf_avg != 0 and non_vaf_avg != 0:
			fc = math.log(ch_vaf_avg/non_vaf_avg,2)
		frow.extend([ch_vaf_avg,non_vaf_avg,fc])
		flines.append(frow)
	for row in vlines:
		vaf_mtx.write('\n' + '\t'.join(map(str,row)))
	vaf_mtx.close()

	for row in flines:
		fc_mtx.write('\n'+'\t'.join(map(str,row)))
	fc_mtx.close()

def parse_array(infile, anno):
	array = {}
	with open(infile) as f:
		lines = f.readlines()
		samples = lines[0].strip().split('\t')[5:]
		samples = map(lambda x: re.sub(r'SCA_',r'SCA-', x), samples)
		for i in range(1, len(lines)):
			line = lines[i].strip()
			col = line.split('\t')
			pos = col[0]; ref = col[1]; alt = col[2]; gene = col[3]
			calls = col[5:]
			var = pos+':'+ref+'>'+alt
			chr = pos.split(':')[0]; p = int(pos.split(':')[1])
			genes = anno[chr]

			for g in genes.keys():
				for po in genes[g]:
					start = po[0]; end = po[1]
					if start <= p and end >= p:
						gene = g
						break
					else: continue

			
			var = chr+':'+str(p)+':'+ref+'>'+alt
		
			if not array.has_key((var,gene)):
				array[(var,gene)] = {}
			for j in range(len(samples)):
				sam = samples[j]
				call = calls[j]
				if call == ref+ref:
					varcheck = 0
				else:
					if ref in call:
						varcheck = 1
					else:
						if '-' in call:
							varcheck = 0
						else:
							if 'I' in call:
								varcheck = 2.5
							elif 'D' in call:
								varcheck = 2.5
							else:
								varcheck = 2
				array[(var,gene)][sam] = [call, varcheck]
	return array				

def parse_chip(infile):
	chips = {}
	with open(infile, 'r') as f:
		lines = f.readlines()
		for i in range(len(lines)):
			line = lines[i].strip()
			col = line.split('\t')
			chr = col[0]; pos = map(int,col[1:3]); gene = col[3].split('.')[0]
			if not chips.has_key(chr):
				chips[chr] = []
			chips[chr].append(pos)
	return chips

sample_file = sys.argv[1]
mtx_file = sys.argv[2]
chip_file = sys.argv[3]
anno_file = sys.argv[4]
outputdir = sys.argv[5]
chip = parse_chip(chip_file)
anno = parse_anno(anno_file)

if not os.path.exists(outputdir): os.makedirs(outputdir)
group = parse_samplefile(sample_file)
mtx = parse_vaf_mtx(mtx_file)

ch_sam = group['CH']; non_ch_sam = group['Non-CH']

array = parse_array(mtx_file, anno)

sam_ids_order = ch_sam+non_ch_sam
lines = [['VAR_ID','CHIP', 'GENE', 'Group1', 'Group2', 'FC']]
lines[0].extend(sam_ids_order)

countdic = {'CH':{}, 'Non-CH':{}}

for var in array.keys():
	samdic = array[var]
	row = [var[0],var[1], 0, 0]
	varid = var[0].split(':')
	chr = varid[0]; pos = int(varid[1])
	chip_pos = chip[chr]
	for cp in chip_pos:
		if pos >= cp[0] and pos <= cp[1]:
			row.insert(1, 'O')
			break
		else: continue
	
	if len(row) != 5:
		row.insert(1,'.')

	for sam in ch_sam:
		call,score = samdic[sam]
		if score > 0:
			row[3] += 1
		row.append(score)
		
		if not countdic['CH'].has_key(sam):
			countdic['CH'][sam] = [0,0,0,0]	
		if score == 0:
			countdic['CH'][sam][0] += 1
		elif score == 1:
			countdic['CH'][sam][1] += 1
		elif score == 2:
			countdic['CH'][sam][2] += 1
		else:
			countdic['CH'][sam][3] += 1
		
	for sam in non_ch_sam:
		call, score = samdic[sam]
		if score > 0:
			row[4] += 1
		row.append(score)

		if not countdic['Non-CH'].has_key(sam):
			countdic['Non-CH'][sam] = [0,0,0,0]
		if score == 0:
			countdic['Non-CH'][sam][0] += 1
		elif score == 1:
			countdic['Non-CH'][sam][1] += 1
		elif score == 2:
			countdic['Non-CH'][sam][2] += 1
		else:
			countdic['Non-CH'][sam][3] += 1

	ch_sum = row[3]; nonch_sum = row[4]
	if ch_sum == 0: fc = 0.0
	elif nonch_sum == 0: fc = 'inf'
	else: fc = (float(ch_sum)/len(ch_sam))/(float(nonch_sum)/len(non_ch_sam))
	row.insert(5,fc)
	lines.append(row)

with open(outputdir + '/array_CDS.CH_vs_NonCH.mtx.txt', 'w') as f:
	for row in lines:
		f.write('\n'+'\t'.join(map(str, row)))

with open(outputdir + '/count_CDS.array_CH_vs_NonCH.txt', 'w') as f:
	f.write('CH\tSamples\tREF\tHET\tHOM')
	for sam in ch_sam:
		counts = countdic['CH'][sam]
		f.write('\nCH_'+sam+'\t'+'\t'.join(map(str,counts)))
	for sam in non_ch_sam:
		counts = countdic['Non-CH'][sam]
		f.write('\nNon-CH_'+sam+'\t'+'\t'.join(map(str,counts)))
