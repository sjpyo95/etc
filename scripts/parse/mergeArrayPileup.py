import sys, os
import vcf

def parse_array(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		samples = header[8:]
		rows = {}
		for line in f:
			line = line.strip().split('\t')
			chrom,pos = line[:2]
			id = line[3]
			af = line[7]
			gene = line[6]
			gt = line[8:]
			rows[(chrom,pos,id,gene,af)] = gt
	return rows, samples

def parse_pileup(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		samples = header[4:]
		rows = {}
		for line in f:
			line = line.strip().split('\t')
			chrom = line[0]; pos = line[1]; chip = line[2]; ref = line[3]
			depths = line[4:]
			rows[(chrom,pos)] = [chip,ref,depths]
	return rows, samples

#def parseVcf(filename, positionlist):
#	reader = vcf.Reader(open(filename, 'r'))
#	ids = {}
#	for record in reader:
#		chrom = record.CHROM; pos = record.POS; id = record.ID
#		if (chrom, pos) in positionlist:
#			ids[(chrom,pos)] = id
#	return ids
		
arrayfile = sys.argv[1]
pileupfile = sys.argv[2]
#dbsnpfile = sys.argv[3]
array, array_samples = parse_array(arrayfile)
pileup, pileup_samples = parse_pileup(pileupfile)
positionlist = array.keys()
positionlist = sorted(positionlist, key = lambda x: (int(x[0]),int(x[1])))
all_samples = array_samples + pileup_samples

with open(sys.argv[3], 'w') as f:
	f.write('CHROM\tPOS\tArrayID\tCHIP\tGENE\tREF\tArrayAF\t'+'\t'.join(all_samples))
	for chrposid in positionlist:
		chrom,pos,id,gene,af = chrposid
		chip,ref,depths = pileup[(chrom,pos)]
		gt = array[chrposid]
		f.write('\n'+'\t'.join([chrom,pos,id,chip,gene,ref,af])+'\t'+'\t'.join(gt) +'\t'+ '\t'.join(depths))

	
	
	
