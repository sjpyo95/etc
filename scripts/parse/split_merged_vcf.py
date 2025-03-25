import sys, os
sys.path.insert(0,'/home/sunyme95/scripts/modules')
import config as cf

def parselist(filename):
	infile = open(filename, 'r')
	lines = infile.readlines(); infile.close()
	samples = []
	for i in range(len(lines)):
		line = lines[i].strip()
		sample = line.split('/')[-1].split('.gvcf')[0]
		samples.append(sample)
	return samples

vcffile = sys.argv[1]
gvcf_list = sys.argv[2]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)
soutputdir = outputdir + '/step2/'
if not os.path.exists(soutputdir): os.makedirs(soutputdir)
samples = parselist(gvcf_list)

for sample in samples:
	outfile = outputdir + '/' + sample +'.vcf.gz'
	cmd = 'bcftools view -s ' + sample + ' ' + vcffile + ' -Oz -o ' + outfile
	cf.run_cmd(cmd)
	cmd2 = 'bcftools norm -m-both -o ' + outputdir + '/' +sample + '.step1.vcf ' + outfile
	cf.run_cmd(cmd2)
	
	cmd3 = 'bcftools norm -f /mnt/mone/PMI/CH/01.Alignment/khyojoo01/001-REFGENOME/genome.fa -o ' + soutputdir + sample + '.step2.vcf ' + outputdir + '/' + sample + '.step1.vcf'
	cf.run_cmd(cmd3)
	cf.run_cmd('rm ' + outputdir + '/' + sample + '.step1.vcf')





