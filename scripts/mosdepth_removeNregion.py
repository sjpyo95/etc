import sys, os
import gzip
def filter_bed_by_fasta(bed_file, fasta_file, gtf_file):
	filtered_lines = []
	with open(fasta_file, 'r') as fasta:
		fasta_data = {}
		current_chrom = ""
		sequence = ""
		for line in fasta:
			line = line.strip()
			if line.startswith(">"):
				if current_chrom:
					fasta_data[current_chrom] = sequence
				current_chrom = line[1:]
				sequence = ""
			else:
				sequence += line
		if current_chrom:
			fasta_data[current_chrom] = sequence
			
	with open(gtf_file, 'r') as gtf:
		gtf_data = {}
		lines = gtf.readlines()
		for line in lines:
			line = line.strip()
			cols = line.strip().split('\t')
			chr = cols[0]; stat = cols[2]; start = int(cols[3]); end = int(cols[4])
			if stat != 'exon': continue
			if not chr in gtf_data:
				gtf_data[chr] = []
			gtf_data[chr].append([start, end])

	with gzip.open(bed_file, 'rt') as bed:
		for line in bed:
			chrom, start, end = line.strip().split()[:3]
			if chrom not in fasta_data:
				continue
			if chrom not in gtf_data:
				continue
			
			sequence = fasta_data[chrom][int(start):int(end)+1]
			if len(set(sequence)) == 1 and 'N' in sequence:
				continue
			if len(filter(lambda x: x[0] <= int(start) or x[1] >= int(end), gtf_data[chrom])) == 0:
				continue
			else:
				filtered_lines.append(line)
	return filtered_lines

def main(bed_file, fasta_file, gtf_file, out_file):
	
#$bed_file = sys.argv[1]
#fasta_file = sys.argv[2]
#out_file = sys.argv[3]
	filtered_lines = filter_bed_by_fasta(bed_file, fasta_file, gtf_file)

	with open(out_file, 'w') as f:
		for line in filtered_lines:
			f.write(line)

if __name__ == "__main__":
	inputdir = sys.argv[1]
	fasta_file = sys.argv[2]
	gtf_file = sys.argv[3]

	samples = os.listdir(inputdir)

	for i in range(len(samples)):
		sample = samples[i]
		print(sample)
		sinputdir = inputdir + sample + '/'
		bed_file = sinputdir + [x for x in os.listdir(sinputdir) if '.regions.bed' in x and '.csi' not in x][0]
		out_file = sinputdir + sample +'.filtered.exon.bed'
		main(bed_file, fasta_file, gtf_file, out_file)
