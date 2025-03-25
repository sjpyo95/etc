import sys, os

# Get the list of VCF files
inputdir = sys.argv[1]
#samples = filter(lambda x: '.vcf' not in x and 'header' not in x, os.listdir(inputdir))
samples = filter(lambda x: '.vcf' in x and not'.idx' in x and not'.stats' in x, os.listdir(inputdir))
#samples = sorted(samples, key=lambda x: int(x.split('-')[1]))
# Create a new VCF file
out_vcf = open(sys.argv[2], "w")

# Write the header line to the new VCF file
out_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+'\t'.join([x.split('.')[0] for x in samples]))

# Loop over the VCF files
newcalls = {}
newinfo = {}
for sample in samples:
#	sinputdir = inputdir + '/'+sample+'/'
#	vcf_file = sinputdir + filter(lambda x: '.vcf' in x and not '.tbi' in x, os.listdir(sinputdir))[0]
	vcf_file = inputdir + '/' + sample
    # Open the VCF file
	vcf = open(vcf_file, "r")
	# Loop over the lines in the VCF file
	for line in vcf:

		# If the line is a header line, skip it
		if line.startswith("#"):
			continue

		# Split the line into fields
		fields = line.split()
		
		# Get the chromosome, position, ID, reference allele, alternate allele, quality score, filter, and INFO fields
		chr = fields[0]
		pos = fields[1]
		id = fields[2]
		ref = fields[3]
		alt = fields[4]
		qual = fields[5]
		filt = fields[6].split(';')
		info = fields[7].split(';')
		call = fields[-1]

		if "1|0" in call or "0|1" in call:
			tmp = call.split(':')
			tmp.pop(-2)
			tmp.pop(-2)
			tmp.pop(-2)
			call = ':'.join(tmp)
		

		varid = ':'.join([chr, pos, ref+'>'+alt])
		if not newcalls.has_key(varid):
			newcalls[varid] = {}
		newcalls[varid][sample] = call

		if not newinfo.has_key(varid):
			newinfo[varid] = [[],[]]
		newinfo[varid][0] += filt; newinfo[varid][1] += info
		newinfo[varid][0] = list(set(newinfo[varid][0])); newinfo[varid][1] = list(set(newinfo[varid][1]))
lines = ''
for varid in newinfo.keys():
	tmp = varid.split(':')
	chr = tmp[0]; pos = tmp[1]; ref = tmp[2].split('>')[0]; alt = tmp[2].split('>')[1]
	samdic = newcalls[varid]
	infos = newinfo[varid]
	filt = infos[0]; info = infos[1]
	line = '\t'.join([chr, pos,'.',ref,alt,'.',';'.join(filt),';'.join(info),'GT:AD:AF:DP:F1R2:F2R1:FAD:SB'])
	calls = 0
	for sample in samples:
		if sample in samdic.keys():
			call = samdic[sample]
			calls+=1
		else:
			call = './.'
			calls+=1
		line += '\t' + call
	lines += '\n' + line

out_vcf.write(lines)

# Close the new VCF file
out_vcf.close()

