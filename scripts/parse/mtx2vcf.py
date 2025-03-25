import sys
#import vcfpy
import csv

mtxfile = sys.argv[1]
outfile = sys.argv[2]


with open(mtxfile, 'rb') as f:
	header = f.readline().strip().split('\t')
	matrix = [line.strip().split('\t') for line in f]
	
#writer = vcfpy.Writer.from_path(outfile, vcfpy.Header(samples=header[5:]))
#
#for row in matrix:
#	chrom_pos, ref, alt = row[0], row[1], row[2]
#	chrom, pos = chrom_pos.split(':')
#	record = vcfpy.Record(
#		CHROM = chrom,
#		POS = int(pos),
#		ID = '.',
#		REF = ref,
#		ALT = [vcfpy.Substitution(alt)],
#		QUAL=None,
#		FILTER=[],
#		INFO={},
#		FORMAT=['GT'],
#		calls=[vcfpy.Call(sample=header[i], data = {'GT':'./.' if row[i] == ref else '0/1' for i in range(5, len(row))}
#		) for i in range(5, len(row))
#		]
#
#	)
#	writer.write_record(record)
#
#writer.close()

with open(outfile, 'w') as f:
	f.write('##fileformat=VCFv4.2\n')
	f.write('''##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">\n##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=3,length=198022430>\n##contig=<ID=4,length=191154276>\n##contig=<ID=5,length=180915260>\n##contig=<ID=6,length=171115067>\n##contig=<ID=7,length=159138663>\n##contig=<ID=8,length=146364022>\n##contig=<ID=9,length=141213431>\n##contig=<ID=10,length=135534747>\n##contig=<ID=11,length=135006516>\n##contig=<ID=12,length=133851895>\n##contig=<ID=13,length=115169878>\n##contig=<ID=14,length=107349540>\n##contig=<ID=15,length=102531392>\n##contig=<ID=16,length=90354753>\n##contig=<ID=17,length=81195210>\n##contig=<ID=18,length=78077248>\n##contig=<ID=19,length=59128983>\n##contig=<ID=20,length=63025520>\n##contig=<ID=21,length=48129895>\n##contig=<ID=22,length=51304566>\n##contig=<ID=X,length=155270560>\n##contig=<ID=Y,length=59373566>\n''')

	f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(header[5:])))

	for row in matrix:
		chrom_pos, ref, alt = row[0], row[1], row[2]
		chrom, pos = chrom_pos.split(':')
		calls = []
		ac = 0
		for i in range(5, len(row)):
			
			if row[i] == ref+ref or row[i] == '.':
				calls.append('0/0')
			elif row[i] == '--':
				calls.append('0/0')
			else:
				if row[i][0] == ref:
					calls.append('1/0')
					ac += 1

				elif row[i][1] == ref:
					calls.append('0/1')
					ac += 1
				else:
					calls.append('1/1')
					ac += 2
		af = float(ac)/(len(row[5:])*2)
		if af >= 0.01:
			af = '{:.3f}'.format(af)
		else:
			af = '{:.3e}'.format(af)
#			print chrom, pos, ref, alt
		alt_value = next((x for x in alt.split(',') if x != '-'), '') if alt != '-' else ''
		if alt_value == '':
			alt_value = ref
		f.write('{}\t{}\t.\t{}\t{}\t.\t.\tAC={};AF={}\tGT\t{}\n'.format(chrom, pos, ref, alt_value, ac, af, '\t'.join(calls)))
