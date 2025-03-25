import sys
import pysam
import vcf

gvcf_file = sys.argv[1]
fasta_file = sys.argv[2]

genome = pysam.FastaFile(fasta_file)

gvcf = vcf.Reader(open(gvcf_file, 'rb'))
sample = gvcf.samples[0]
coverage = {}

for record in gvcf:
	chrom = record.CHROM
	pos = record.POS
	ref_base = record.REF
	all_bases = record.ALT
	data = record.genotype(sample).data
	depth = data.DP
	
#	if 'N' in genome.fetch(chrom, pos-1, pos):
#		continue
	
	## If no alternate alleles
	if len(all_bases) == 1:
		coverage[(chrom, pos)] = depth

	## Calculate the coverage for positions with alternate alleles
	else:
		alt_bases = all_bases[:-1]
		alt_ads = data.AD[:-2]
		ref_ad = data.AD[-2]
		
		if len(alt_bases) == 1:
			alt_base = alt_bases[0]
			alt_ad = alt_ads[0]
			if len(ref_base) == len(alt_base):
				for p in range(len(ref_base)):		
					coverage[(chrom, pos+p)] = alt_ad + ref_ad
			
			# Deletion
			elif len(ref_base) > len(alt_base):
				for p in range(len(ref_base)):
					if p < len(alt_base):
#						print(chrom, pos+p, ref_base, alt_base, alt_ad, ref_ad, alt_ad+ref_ad, depth)
						coverage[(chrom,pos+p)] = alt_ad + ref_ad
					else:
#						print(chrom, pos+p, ref_base, alt_base, alt_ad, ref_ad, ref_ad, depth)
						coverage[(chrom, pos+p)] = ref_ad
		## Insertion1
			elif len(ref_base) < len(alt_base):
				for p in range(len(alt_base)):
					if p < len(ref_base):
#						print(chrom, pos+p, ref_base, alt_base, alt_ad, ref_ad, alt_ad+ref_ad)
						coverage[(chrom, pos+p)] = alt_ad+ref_ad
					else:
#						print(chrom, pos+p, ref_base, alt_base, alt_ad, ref_ad, alt_ad)
						coverage[(chrom, pos+p)] = alt_ad
#			else:
#				print(chrom, pos, ref_base, alt_base, alt_ad, ref_ad, alt_ad+ref_ad, depth)
		else:
			for i in range(len(alt_bases)):
				if (chrom, pos) in coverage: continue
				alt_base = alt_bases[i]
				alt_ad = alt_ads[i]
				ref_len = len(ref_base)
				alt_len = len(ref_base)
				for p in range(max(ref_len, alt_len)):
					if (chrom, pos+p) not in coverage:
						coverage[(chrom, pos+p)] = 0

					if ref_len > p:
						coverage[(chrom, pos+p)] += ref_ad
					
					if alt_len > p:
						
						coverage[(chrom, pos+p)] += alt_ad

						
#				print(ref_base, alt_base,alt_ad, data.AD)
print(coverage)
print ('Chromosome\tPosition\tRef Base\t Alt Base\tVarStat\tAD\tDP')
for key, value in coverage.items():
	print('\t'.join(map(str,key)) + '\t' + '\t'.join(map(str,value)))

