import sys, os

def parse_mpileup_line(line):
	line = line.strip().split("\t")
	chrom = line[0]
	pos = line[1]
	ref = line[2]
	depth = int(line[3])
	bases = line[4].upper()
	counts = {}
	i = 0
	while i < len(bases):
		if bases[i] == '^':
			i += 2
		elif bases[i] == '*':
			i += 1
		elif bases[i] == '$':
			i += 1
		elif bases[i] == '+':
			j = i + 1
			while j < len(bases) and bases[j].isdigit():
				j += 1
			indel_len = int(bases[i+1:j])
			indel_bases = '+'+str(indel_len)+bases[j:j+indel_len]
			counts[indel_bases] = counts.get(indel_bases, 0) + 1
			i = j+indel_len
		elif bases[i] == '-':
			j = i + 1
			while j < len(bases) and bases[j].isdigit():
				j += 1
			indel_len = int(bases[i+1:j])
			indel_bases = '-'+str(indel_len)+bases[j:j+indel_len]
#			indel_bases = ref[int(pos)+i-1:int(pos)+i-1+indel_len]
			counts[indel_bases] = counts.get(indel_bases,0)+ 1
			i = j + indel_len
		else:
			base = bases[i]
			if base == '.' or base == ',' : base = 'REF('+ref+')'
			
			counts[base] = counts.get(base, 0) + 1
			i += 1
	
	return chrom, pos,ref,counts,depth

def main(inputfile, outputfile):
	outfile = open(outputfile, 'w')
	outfile.write('CHROM\tPOS\tREF\tDEPTH')
	with open(inputfile, 'r') as f:
		for line in f:
			chrom, pos, ref, counts, depth = parse_mpileup_line(line)
#			print line
			counts = sorted(counts.items(), key=lambda x:x[1], reverse = True)
			maf = 0
			if len(counts) > 1:
				mjcount = float(counts[0][1]); mncount = float(counts[1][1])
				maf = mncount/(mjcount+mncount)
			counts.append(('MAF',maf))
#			if 'REF' not in counts[0][0]: print counts
#			print('{}\t{}\t{}\t{}'.format(chrom, pos, ref, ';'.join(['{}:{}'.format(k,v) for k,v in counts])))
			outfile.write('\n{}\t{}\t{}\t{}'.format(chrom, pos, ref, ';'.join(['{}:{}'.format(k,v) for k,v in counts])))
	outfile.close()

if __name__ == '__main__':
	py, inputdir, outputdir = sys.argv
	if not os.path.exists(outputdir): os.makedirs(outputdir+'/')
	files = [inputdir +'/'+ x for x in filter(lambda x: '.pilup' in x, os.listdir(inputdir))]
	for i in range(len(files)):
		infile = files[i]
		sampleName = infile.split('/')[-1].split('.pilup')[0]
		outfile = outputdir +'/'+ sampleName + '.depth.txt'
		main(infile, outfile)
		
	
	
