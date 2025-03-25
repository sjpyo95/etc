import sys, os

def parse_pileup(filename):
	with open(filename, 'r') as f:
		dic = {}
		lines = f.readlines()
		for i in range(1,len(lines)):
			line = lines[i].strip()
			col = line.split('\t')
			chrom = col[0]; pos = col[1]
			ref = col[2]; alt = col[3]
			depth = col[4]
			dic[(chrom, pos)] = float(depth)

	return dic		

def compare_depth(pile1, pile2, outfile):
	with open(outfile, 'w') as f:
		for pos in pile1.keys():
			if pos in pile2.keys():
				depth1 = pile1[pos]
				depth2 = pile2[pos]
				dif = depth1 - depth2
				f.write('\t'.join(pos)+'\t'+'\t'.join(map(str, [depth1, depth2, dif]))+'\n')			



def main(pile1file, pile2file, outfile):
	pile1 = parse_pileup(pile1file)
	pile2 = parse_pileup(pile2file)

	compare_depth(pile1, pile2, outfile)



if __name__ == "__main__":
	inputdir1 = sys.argv[1]
	inputdir2 = sys.argv[2]
	ill = filter(lambda x: '.pileup' in x, os.listdir(inputdir1))
	ult = filter(lambda x: '.pileup' in x, os.listdir(inputdir2))

	for i in ill:
		for u in ult:
			if u.split('.')[0] in i:
				outfile = outputdir + u.split('.')[0] + '.Depth_Difference.txt'
				compare_depth(inputdir1+i, inputdir2+u,outfile)
				break

