import sys, os

with open(sys.argv[1], 'r') as f:
	lines = f.readlines()

with open(sys.argv[2], 'w') as f:
	f.write('##fileformat=VCFv4.2\n')
	f.write('##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=3,length=198022430>\n##contig=<ID=4,length=191154276>\n##contig=<ID=5,length=180915260>\n##contig=<ID=6,length=171115067>\n##contig=<ID=7,length=159138663>\n##contig=<ID=8,length=146364022>\n##contig=<ID=9,length=141213431>\n##contig=<ID=10,length=135534747>\n##contig=<ID=11,length=135006516>\n##contig=<ID=12,length=133851895>\n##contig=<ID=13,length=115169878>\n##contig=<ID=14,length=107349540>\n##contig=<ID=15,length=102531392>\n##contig=<ID=16,length=90354753>\n##contig=<ID=17,length=81195210>\n##contig=<ID=18,length=78077248>\n##contig=<ID=19,length=59128983>\n##contig=<ID=20,length=63025520>\n##contig=<ID=21,length=48129895>\n##contig=<ID=22,length=51304566>\n##contig=<ID=X,length=155270560>\n##contig=<ID=Y,length=59373566>\n')
	f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n')
	f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
	for i in range(len(lines)):
		line = lines[i].strip().split('\t')
		tmp = line[0].split('_')
		af = float(line[1])
		chr = tmp[0].split('chr')[1]; pos = tmp[1]; ref = tmp[2]; alt = tmp[3]
		f.write('\t'.join([chr,pos,'.',ref,alt,'.','.','AF='+str(af)])+'\n')



