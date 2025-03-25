import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf

inputdir = sys.argv[1]
#vcffile = sys.argv[1]
outputdir = sys.argv[2]
samples = filter(lambda x: '.vcf' not in x and '.txt' not in x and 'logs' not in x, os.listdir(inputdir))

#soutputdir = outputdir + '/split/'
#if not os.path.exists(soutputdir):os.makedirs(soutputdir)
#cf.run_cmd('bcftools +split ' + vcffile + ' -o ' + soutputdir)


#inputdir = soutputdir


#files = filter(lambda x: 'srt.vcf' in x, os.listdir(inputdir))
#for filename in os.listdir(inputdir):
for sample in samples:
	sinputdir = inputdir + '/' + sample + '/'
	samfile = sinputdir + filter(lambda x: '.vcf' in x and not '.srt' in x, os.listdir(sinputdir))[0]
	srtfile = sinputdir + sample + '.srt.vcf'
	cf.run_cmd('grep -v "0.0:" ' + samfile + ' > ' + srtfile)
#	if filename.endswith('.vcf'):
#		samfile = os.path.join(inputdir, filename)
#		srtfile = os.path.join(inputdir, filename.replace('.vcf', '.srt.vcf'))
#		cf.run_cmd('grep -v "0.0:" ' + samfile + ' > ' + srtfile)

ssoutputdir = outputdir + '/mafs/'
if not os.path.exists(ssoutputdir):os.makedirs(ssoutputdir)

#for filename in os.listdir(inputdir):
for sample in samples:
	sinputdir = inputdir + '/' + sample + '/'
	srtfile = sinputdir + sample + '.srt.vcf'
	outfile = ssoutputdir + sample + '.maf'
	cf.run_cmd('vcf2maf.pl --input-vcf ' + srtfile + ' --output-maf ' + outfile + ' --ref-fasta /mnt/mone/PMI/CH/Reference/Genome/bwa-mem/genome.fa --inhibit-vep --tumor-id ' + sample)
#	if filename.endswith('srt.vcf'):
#		srtfile = os.path.join(inputdir, filename)
#		outfile = os.path.join(ssoutputdir, filename.replace('.srt.vcf', '.maf'))
#		id = filename.split('.srt.vcf')[0]
#		cf.run_cmd('vcf2maf.pl --input-vcf ' + srtfile + ' --output-maf ' + outfile + ' --ref-fasta /mnt/mone/PMI/CH/Reference/Genome/bwa-mem/genome.fa --inhibit-vep --tumor-id ' + id)

cf.run_cmd('head -n 2 ' + outfile + ' > ' + outputdir + '/head.tmp')

cf.run_cmd('cat ' + ssoutputdir + '*.maf | grep -v "#" | grep -v "Hugo_Symbol" > ' + outputdir + '/merged.maf.tmp')
cf.run_cmd('''awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k5,5V -k6,6n"}' ''' + outputdir +'/merged.maf.tmp > ' + outputdir + '/merged.maf.tmp2')
cf.run_cmd('cat ' + outputdir + '/head.tmp ' + outputdir + '/merged.maf.tmp2 > ' + outputdir + '/'+ 'merged.maf')

cf.run_cmd('rm ' + outputdir + '/*.tmp*')
#cf.run_cmd('rm ' + ssoutputdir)
#cf.run_cmd('rm ' + soutputdir)
