import sys, os

def parseMtx(filename):
	with open(filename,'r') as f:
		lines = f.readlines()
		col = lines[0]
		mtx = dict()
		for i in range(1, len(lines)):
			line = lines[i].strip()
			cols = line.split('\t')
			varid = cols[0]
			pos = varid.split(':')[0,1]
			mtx[pos] = line
			
	return col, mtx

def parseCHIP(filename):
	with open(filename, 'r') as f:
		lines = f.readlines()
		genes = dict()
		for i in range(len(lines)):
			line = lines[i].strip()
			cols = line.split('\t')
			chr = cols[0]; start = cols[1]; end = cols[2]
			id = [chr, start, end]
			gene = cols[3].split('.')[0]
			if not genes.has_key(gene):
				genes[gene] = []
			genes[gene].append(id)
	return genes
	
mtxfile = sys.argv[1]
col, mtx = parseMtx(mtxfile)

chipgenefile = sys.argv[2]
chipgenes = parseCHIP(chipgenefile)

for gene in chipgenes.keys():
	positions = chipgenes[gene]
	chr = positions[0][0]
	start = min([int(x[1]) for x in positions])
	end = max([int(x[2]) for x in positions])
	filter(lambda x:
