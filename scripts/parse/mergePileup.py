import sys, os
import re
def parsePileup(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		pileup = {}
		for line in f:
			chrom,pos,ref,depth = line.strip().split('\t')
			depth_info = depth.split(';')
			refdepth = int(depth_info[0].split(':')[1])
			alts = depth_info[1:-1]
			altdepths = [re.split(':', x) for x in alts]
			if len(alts) < 1:
				alts = [ref]
#			maxmaf = float(depth_info[-1].split(':')[1])*100
			for x in altdepths:
				x[1] = float(x[1])/refdepth
				if x[1] < 0.3:
					altdepths.remove(x)
			pileup[(chrom,pos)] = [ref,refdepth, alts, altdepths]

	return pileup

def parsePOS(filename):
	with open(filename, 'r') as f:
		header = f.readline().strip().split('\t')
		chip = []
		for line in f:
			line = line.strip().split('\t')
			chrom = line[0]; pos = line[1]
			chip.append((chrom,pos))
	return chip
	

def main():
	inputdir = sys.argv[1]
	filenames = filter(lambda x: 'depth.txt' in x , os.listdir(inputdir))
#	interPos = parsePOS(sys.argv[2])
	chipPOS = parsePOS(sys.argv[2])
	all_pileup = {}
	samples = []
	for i in range(len(filenames)):
		depthfile = inputdir + '/' + filenames[i]
		sample = re.sub(r'IDT_',r'',filenames[i])
		sample = re.sub(r'.depth.txt',r'',sample)
		samples.append(sample)

	samples = sorted(samples, key = lambda x: int(x.split('-')[1]))
	pileup = parsePileup(depthfile)


	with open(sys.argv[3], 'w') as f:
		f.write('VARID\tCHIP\tREF\tALT\t'+'\t'.join(samples))

		for p in all_pileup.keys():
			chrom, pos, ref = p
			if (chrom,pos) in chipPOS:
				chip = 'O'
			else: 
				chip = 'X'
			f.write('\n'+'\t'.join([chrom+':'+pos,chip,ref]))
			for sample in samples:
				depth = all_pileup[p][sample]
				f.write('\t'+depth)
			

if __name__=='__main__':
	main()
###Command: python2.7 mergePileup.py depthdir/ intersect.CHROM_POS.txt known_CHIPs.CHROM_POS.txt outfile.txt
