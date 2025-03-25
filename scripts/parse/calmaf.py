import sys, os


with open(sys.argv[1], 'r') as f:
	header = f.readline().strip().split('\t')
	samples = header[4:]
#	print '\t'.join(header)
	print 'Sample\tCHROM\tPOS\tCHIP\tMAF\tMAJ_BASE\tMIN_BASE\tMAJ_COUNT\tMIN_COUNT_SUM'
	for line in f:
		newline = ''
		line = line.strip().split('\t')
		chrom = line[0]; pos = line[1]; chip = line[2]; ref = line[3]
		depths = line[4:]
		newline += '\t'.join([chrom,pos,chip,ref])
		
		for i in range(len(samples)):
			sample = samples[i]
			depth = depths[i]
			if len(depth.split(';')) > 1:
				major = float(depth.split(';')[0].split(':')[1])
				major_allele = depth.split(';')[0].split(':')[0]
#				minor = sum(float(x.split(':')[1]) for x in depth.split(';')[1:])
				minor = float(depth.split(';')[1].split(':')[1])
#				minor_alleles = [x.split(':')[0] for x in depth.split(';')[1:]]
				minor_alleles = depth.split(';')[1].split(':')[0]
				maf = round(minor/(major+minor), 6)
				if maf > 0.4: continue 
				if  maf < 0.05: continue
				print '\t'.join([sample,chrom,pos,chip,str(maf),major_allele,minor_alleles, str(major),str(minor)])

			else:
				maf = 0
#				print '\t'.join([sample,chrom,pos,chip,str(maf),major_allele, '.', str(major), '.'])
			newline+= '\t'+depth+';MAF:'+str(maf)
#		print newline
#			depth += ';MAF:'+str(maf)


