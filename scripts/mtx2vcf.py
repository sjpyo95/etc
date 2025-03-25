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
	f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(header[5:])))

	for row in matrix:
		chrom_pos, ref, alt = row[0], row[1], row[2]
		chrom, pos = chrom_pos.split(':')
		calls = []
		for i in range(5, len(row)):
			if row[i] == ref:
				calls.append('./.')
			elif row[i] == '-':
				calls.append('.')
			else:
				calls.append('0/1')
		alt_value = next((x for x in alt.split(',') if x != '-'), '') if alt != '-' else ''
		if alt_value == '':
			alt_value = ref
		f.write('{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t{}\n'.format(chrom, pos, ref, alt_value, '\t'.join(calls)))
