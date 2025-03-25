import sys
import vcf
import pysam

gvcf_file = sys.argv[1]
fasta_file = sys.argv[2]

vcf_reader = vcf.Reader(open(gvcf_file,'rb'))

genome = pysam.FastaFile(fasta_file)

#total_pos = 0
#cover_pos = 0

#chr_total_pos = {}
#chr_cover_pos = {}
#chr_n_base_len = {}
uncovered_length = 0
chr_uncovered_len = {}
for chrom in vcf_reader.contigs:
	chrom_len = genome.get_reference_length(chrom)
	chr_uncovered_len[chrom] = chrom_len
	start = 0

	for record in vcf_reader.fetch(chrom):
		end = record.POS - 1
		covered_length = end - start + 1
		chr_uncovered_len[chrom] -= covered_length
		uncovered_length += covered_length

		start = record.POS

	uncovered_length += chrom_len - start

total_length = sum(chr_uncovered_len.values())
percent_uncovered = (uncovered_length/total_length)*100
chr_per_uncovered = {chrom: (uncovered_lengths/chrom_len)*100 for chrom, uncovered_lengths in chr_uncovered_len.items()}

for chrom in chr_per_uncovered.keys():

	print(f'{chrom}\t{chr_uncovered_len[chrom]}\t{chr_uncovered_len[chrom]}\t{chr_per_uncovered[chrom]:.2f}%')
print(f'Total\t{total_length}\t{uncovered_length}\t{percent_uncovered:.2f}%')

#	curr_n_base_len = 0
#	for pos in range(chrom_len):
#		ref_base = genome.fetch(chrom, pos+1, pos+2).upper()
#		if ref_base == 'N':
#			curr_n_base_len += 1
#		else:
#			if curr_n_base_len > 0:
#				chr_n_base_len[chrom] = chr_n_base_len.get(chrom,0) + curr_n_base_len
#				curr_n_base_len = 0
#				
#		total_pos += 1
#	
#
#	chr_total_pos[chrom] = chrom_len - chr_n_base_len.get(chrom, 0)
#	chr_cover_pos[chrom] = 0
#
#	for record in vcf_reader.fetch(chrom):
#		pos = record.POS -1
#		ref_base = genome.fetch(chrom, pos+1,pos+2).upper()
#
#		if ref_base != 'N':
#			alt = record.ALT
#			if len(alt) == 1:
#				cover_pos += 1
#				chr_cover_pos[chrom] += 1
#			else:
#				alt_len = [len(x) for x in alt[:-1]]

#percent_covered = (float(cover_pos)/total_pos)* 100
#
#chr_per_cov = {}
#for chrom in chr_total_pos.keys():
#	chr_per_cov[chrom] = (float(chr_cover_pos[chrom])/chr_total_pos[chrom])*100
#
#print('Chr\tTotal\tCovered\tPercentage')
#for chrom in chr_per_cov.keys():
#	print('{}\t{}\t{}\t{}'.format(chrom, chr_total_pos[chrom], chr_cover_pos[chrom], chr_per_cov[chrom]))
#
#print('Total\t{}\t{}\t{}'.format(total_pos, cover_pos, percent_covered))

