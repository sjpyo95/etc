import vcf
from pyfaidx import Fasta
import sys
import gzip

vcf_file = sys.argv[1]
fasta_file = sys.argv[2]

vcfread = vcf.Reader(gzip.open(vcf_file, 'r'))
fasta = Fasta(fasta_file)
depths = {}

for chr in vcfread.contigs:
	depths[chr] = {}

	for rec in vcfread.fetch(chr):
		ref_base = fasta[chr][rec.POS-1]
		
		if ref_base == 'N': continue

		if record.is_snp:
			depth = rec.INFO['DP']
		elif rec.is_indel:
			if len(rec.REF) < len(rec.ALT[0]):
				depth = rec.INFO['DP']
				for i in range(len(rec.ALT[0]) - len(rec.REF)):
					depths[chr][rec.POS + i] = depth
			else:
				depth = rec.INFO['DP']
				for i in range(len(rec.REF) - len(rec.ALT[0])):
					depths[chr][rec.POS] = depth
					rec.POS += 1
		else:
			depth = 0


		depths[chr][rec.POS] = depth

for chr in depths.keys():
	print (chr,'\t',)
	poss = depths[chr]
	for pos in poss.keys():
		depth = poss[pos]
		print (str(pos),'\t',str(depth))

		
	
