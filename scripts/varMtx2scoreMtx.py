import csv
import sys

filename = sys.argv[1]
outfilename = sys.argv[2]
samstartpos = int(sys.argv[3])
with open(filename, 'rb') as f:
	reader = csv.reader(f, delimiter='\t')
	header = next(reader)
	data = []
	for row in reader:
		data.append(row)

for i in range(len(data)):
	for j in range(samstartpos, len(data[i])):
		if data[i][j] == '.':
			data[i][j] = '0'
		elif data[i][j][0] == data[i][1]:
			data[i][j] = '2'
		elif data[i][j][1] == data[i][1]:
			data[i][j] = '2'
		else:
			data[i][j] = '4'
outfile = open(outfilename, 'w')
outfile.write('\t'.join(header))
for row in data:
	outfile.write('\n' + '\t'.join(row))
outfile.close()
