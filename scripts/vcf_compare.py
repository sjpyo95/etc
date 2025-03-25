import vcf
import pysam
import sys

vcffile1 = sys.argv[1]
vcffile2 = sys.argv[2]

vcf1 = vcf.Reader(open(vcffile1,'rb'))
vcf2 = vcf.Reader(open(vcffile2,'rb'))
#genome = pysam.FastaFile(fastafile)

vcf1dic = {}
for record in vcf1:
	chrom = record.CHROM; pos = record.POS
	refbase = record.REF
	altbases = record.ALT
	if None in altbases:
		dp = record.samples[0]['DP']
	else:
		dp = sum(record.samples[0]['AD'])
	vcf1dic[(chrom, pos, refbase)] = dp


vcf2dic = {}
for record in vcf2:
	chrom = record.CHROM; pos = record.POS
	refbase = record.REF
	altbases = record.ALT
	dp = record.samples[0]['DP']
	vcf2dic[(chrom, pos, refbase)] = dp

common_pos = list(set(vcf1_pos) & set(vcf2_pos))
vcf1_only = list(set(vcf1_pos) - set(vcf2_pos))
vcf2_only = list(set(vcf2_pos) - set(vcf1_pos))


for pos in common_pos:
	print('{}\t{}\t{}\t{}\t{}'.format(pos[0], pos[1], pos[2], vcf1[pos], vcf2[pos]))

print('\n')
for pos in vcf1_only:
	print('{}\t{}\t{}\t{}\t{}'.format(pos[0], pos[1], pos[2], vcf1[pos], 0))

print('\n')
for pos in vcf2_only:
	print('{}\t{}\t{}\t{}\t{}'.format(pos[0], pos[1], pos[2], 0, vcf2[pos]))


