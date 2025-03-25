import sys
sys.path.insert(0,'/home/sunyme95/scripts/modules/')
import os
import config as cf

fasta = '/home/sunyme95/Reference/FASTA/hg38.fa'
#bedfile = '/mnt/mone/PMI/CH/Reference/On_Target_Bed/Targets-XGEN.66F814ECD1264A978BE2D43A40BDAB94.bed'

#script, inputdir, mapq, outputdir = sys.argv

script, inputdir, outputdir = sys.argv

if not os.path.exists(outputdir): os.makedirs(outputdir)

bamfiles = [inputdir + x for x in filter(lambda x: '.bam' in x and '.bai' not in x, os.listdir(inputdir))]

for i in range(len(bamfiles)):
	bamfile = bamfiles[i]
	samplename = bamfile.split('/')[-1].split('_S')[0]
	outfile = outputdir + samplename+'.pileup'
	queue = cf.queue(2)
	samtools = 'samtools mpileup ' + bamfile + ' -f ' + fasta + ' -o ' + outfile
#	samtools = 'samtools mpileup ' + bamfile + ' -q ' + mapq + ' -f ' + fasta + ' --positions ' + bedfile + ' -o ' + outfile
#	cf.run_cmd(samtools)
	cf.qsub_execute(samplename+'_pileup', queue, samtools, 14, outputdir)
#	cf.qsub_time(20, 5)
