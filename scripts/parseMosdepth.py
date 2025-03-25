import sys, os
import gzip

def parseDepth(infile):
	newlines = []
	with gzip.open(infile, 'rb') as f:
		for line in f:
			line = line.decode('utf-8').strip()
			chr, start, end, depth = line.split('\t')
			name = chr+':'+start+'-'+end
			newlines.append(line+'\t'+name)
	return newlines

def write(lines, outfile):
	with open(outfile, 'w') as f:
		f.write('\n'.join(lines))

def main():
	inputdir = sys.argv[1]
	outputdir=  sys.argv[2]
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	samples = [x for x in os.listdir(inputdir) if '-' in x]

	for i in range(len(samples)):
		sample = samples[i]
		sinputdir = inputdir + sample + '/'
		infile = sinputdir + [x for x in os.listdir(sinputdir) if '.regions.bed.gz' in x][0]
		newlines = parseDepth(infile)
		outfile = outputdir + sample + '.new.regions.bed'
		print(outfile)
		write(newlines, outfile)
	
if __name__=="__main__":
	main()
