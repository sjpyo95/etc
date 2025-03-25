import argparse
import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')

import config as cf
import time
import commands
ref = "/mnt/mone/PMI/CH/Reference/Genome/genome.fa"

def groupFastq(inputdir, sep):
#	f = ['105', '387', '338', '336', '275', '284', '295', '309', '373']
#	303

	fastqs = [inputdir + x for x in filter(lambda x: x in f and not '.txt' in x and '.fastq' in x or '.fq' in x, os.listdir(inputdir))]
	samples = list(set([x.split(inputdir)[-1].split(sep)[0] for x in fastqs]))
	fastqss = []
	for sample in samples:
		sam_fastqs = [s for s in fastqs if sample in s]
		sam_fastqs.sort

		fastqss.append(sam_fastqs)
	return samples, fastqss

### FASTQC ###
def run_fastqc(fastqs, outputdir, sample, queue, tN, trim):
	if trim:
		soutputdir = outputdir + '/trim-fastqc/'
	else:
		soutputdir = outputdir + '/fastqc/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	for fastq in fastqs:
		fastqc = 'fastqc -t ' + str(tN) + ' -o ' + soutputdir + ' ' + fastq
		cf.qsub_execute('fastqc.'+sample, queue, fastqc, tN, outputdir)

	return soutputdir

### Trim-galore ###
def run_trimGalore(fastqs, outputdir, sample, queue, tN):
	soutputdir = outputdir + '/trim-galore/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
#	adapt = '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
	clip = ''
	trimgalore = 'trim_galore --paired -j '+ tN + ' ' + clip + ' -o ' + soutputdir + ' ' + fastqs[0] + ' ' + fastqs[1]
	cf.qsub_execute('trim-galore.'+sample, queue, trimgalore, tN, outputdir)

	return soutputdir

### BWA-MEM ###
def run_bwaMem(fastqs, sample, outputdir, queue, tN):
	soutputdir = outputdir + '/sam/' + sample + '/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	bwamem = 'bwa mem -M -t ' + tN + ' ' + ref + ' ' + fastqs[0] + ' ' + fastqs[1] + ' > ' + soutputdir + '/'+ sample + '.sam'
	cf.qsub_execute('bwa-mem.'+sample, queue, bwamem, tN, outputdir)

### Sort SAM to BAM ###
def run_sortSam(samfile, sample, outputdir, queue, tN):
#	soutputdir1 = outputdir+ '/bam/'
#	if not os.path.exists(soutputdir1): os.makedirs(soutputdir1)
	sortSam = 'python2.7 /home/sunyme95/scripts/fq-preprocess/sortSam.py ' + samfile +  ' ' + outputdir + ' ' + sample + ' ' + tN
	cf.qsub_execute('sortSam.'+sample, queue, sortSam, tN, outputdir)

def split2Chr(bamfile, outputdir):
	"""Separates a BAM file into chromosome-specific BAM files using samtools."""
	commands.getstatusoutput("samtools index" + bam_file)  # Ensure BAM file is indexed
	sample = bam_file.split('.')[0]
	with open(bam_file, "rb") as f:
		headers = f.readlines(10000)  # Read header lines
		for header in headers:
			if header.startswith("@SQ"):
				chrom = header.split()[1].replace("SN:", "")  # Extract chromosome name
				output_file = outputdir + sample + '_' + chrom + ".bam"
				commands.getstatusoutput("samtools view -b" + bam_file + ' ' +  chrom + ' > ' + output_file)


#	sam2bam = 'samtools view -@ ' + tN + ' -bS ' + samfile + ' | samtools sort -@ ' + tN + ' -n -o ' + soutputdir1 + sample + '.bam -O BAM'
#	cf.qsub_execute('sam2bam.'+sample, queue, sam2bam, tN, outputdir)
#	cf.qsub_time(1,1)
#
#	bamfix = 'samtools fixmate -@ ' + tN + ' -O BAM -m ' + soutputdir1 + sample + '.bam ' + soutputdir1 + sample + '.fix.bam'
#	cf.qsub_execute('bamfix.'+sample, queue, bamfix, tN, outputdir)
#	cf.qsub_time(1,1)
#
#	soutputdir2 = outputdir + '/sortBam/'
#	if not os.path.exists(soutputdir2):os.makedirs(soutputdir2)
#	fixsort = 'samtools sort -@ ' + tN + ' -o ' + soutputdir2 + sample + '.srt.bam ' + soutputdir1 + sample + '.fix.bam'
#	cf.qsub_execute('fixsort.'+sample, queue, fixsort, tN, outputdir)
#	cf.qsub_time(1,1)
#
#	soutputdir3 = outputdir + '/dupMark/'
#	if not os.path.exists(soutputdir3): os.makedirs(soutputdir3)
#
#	dupmark = 'samtools markdup -@ ' + tN + ' ' + soutputdir2 + sample + '.srt.bam ' + soutputdir3 + sample + '.srt.mrk.bam'
#	cf.qsub_execute('dupmark.'+sample, queue, dupmark, tN, outputdir)
#	cf.qsub_time(1,1)
#
#	index = 'samtools index -@ ' + tN + ' ' + soutputdir3 + sample + '.srt.mrk.bam'
#	cf.qsub_execute('index.'+sample, queue, index, tN, outputdir)
	
if __name__ == '__main__':
	prs = argparse.ArgumentParser(description = "Run pipeline of pre-processing DNA-seq data fastq to BAM in multi-samples")
	prs.add_argument('-I', '--inputfile', help = 'Input directory of fastq files', dest = 'inputdir')
	prs.add_argument('-S', '--samplename', help = 'specific sample name', dest = 'sample')
	prs.add_argument('-O', '--outputdir', help = 'Output directory for BAM files', dest = 'outputdir')
	prs.add_argument('-p', '--process', help = 'Specific process to begin. if start to beginning: all', dest = 'process')
	prs.add_argument('-t', '--thread', help = 'Thread number for qsub', default = 4, dest = 'tN')
	prs.add_argument('-j', '--job', help = 'Job number of multi run', default = 4, type = int, dest = 'jN')
	prs.add_argument('-q', '--queue', help = 'qsub number 2~8, 0', default = 0, type = int, dest = 'qN')
	prs.add_argument('-s', '--separate', help = 'standard string to separate samples (e.g., _R)', dest = 'sep')
	
	ag = prs.parse_args()
	inputdir = ag.inputdir; outputdir = ag.outputdir; process = ag.process; tN = ag.tN; jN = ag.jN; qN = ag.qN; sep = ag.sep
	sample = ag.sample
	queue = cf.queue(qN)

	if not inputdir.endswith('/') : inputdir = inputdir + '/'
	if not outputdir.endswith('/'): outputdir = outputdir + '/'

#	samples, fastqss = groupFastq(inputdir, sep)
	fastqs = filter(lambda x: sample in x and '.fastq' in x or '.fq' in x, os.listdir(inputdir))
	fastqs.sort
	
	all_process = ['all','fastqc','trim-galore','trim-fastqc','bwa-mem', 'sortSam']
	stpN = all_process.index(process)
	all_process = all_process[stpN:]


	if 'fastqc' in all_process:
		print '\n**FASTQC**\n'
#		for i in range(len(fastqss)):
#			fastqs = fastqss[i]
#			sample = filter(lambda x: x in fastqs[0] and '.txt' not in x, samples)[0]
		fastqc_outdir = run_fastqc(fastqs, inputdir, sample, queue,tN, False)
		cf.qsub_time(jN, 1, 'fastqc')
		cf.qsub_time(1,1, 'fastqc')

	if 'trim-galore' in all_process:
		print '\n**Trim-galore**\n'
#		for i in range(len(fastqss)):
#			fastqs = fastqss[i]
#			sample = filter(lambda x: x in fastqs[0] and '.txt' not in x, samples)[0]
		trim_outdir = run_trimGalore(fastqs, inputdir, sample, queue, tN)
		#cf.qsub_time(jN, 1, 'trim')
		cf.qsub_time(1,1, 'trim')
	fastqs = filter(lambda x: sample in x and '.fastq' in x or '.fq' in x, os.listdir(inputdir+'/trim-galore/'))
	fastqs.sort
#	samples, trim_fastqss = groupFastq(inputdir+'/trim-galore/', sep)
	if 'trim-fastqc' in all_process:
		print '\n**Trim-FASTQC**\n'
#		for i in range(len(trim_fastqss)):
#			fastqs = trim_fastqss[i]
#			sample = filter(lambda x: x in fastqs[0] and '.txt' not in x, samples)[0]
		trmFastqc_outdir = run_fastqc(fastqs, inputdir, sample, queue, tN, True)
		cf.qsub_time(jN,1, 'fastqc')
		cf.qsub_time(1,1, 'fastqc')

	if 'bwa-mem' in all_process:
		print '\n**BWA-MEM**\n'
#		for i in range(len(trim_fastqss)):
#			fastqs = trim_fastqss[i]
#			sample = filter(lambda x: x in fastqs[0] and '.txt' not in x, samples)[0]
		run_bwaMem(fastqs, sample, outputdir, queue, tN)
		cf.qsub_time(jN,1, 'bwa')
		cf.qsub_time(1,1, 'bwa-me')
	
	sam_dir = outputdir + '/sam/' + sample + '/'
	samfile = sam_dir + '*.sam'
	if 'sortSam' in all_process:
		print '\n**Sort SAM**\n'
#		samfiles = [sam_dir + x for x in filter(lambda x: '.sam' in x, os.listdir(sam_dir))]
#		for samfile in samfiles:
#			sample = filter(lambda x: x in samfile, samples)[0]
		run_sortSam(samfile, sample, outputdir, queue, tN)
		cf.qsub_time(jN,1,'sortS')
		cf.qsub_time(1,1, 'sortS')

	dpmrkdir = outputdir + '/dpmrk/' + sample + '/'
	bamfile = dpmrkdir + '*.bam'
	if 'splitChr' in all_process:
		print '\n**Split to Chromosomes**\n'
		splitdir = dpmrkdir + '/chrs/'
		if not os.path.exists(splitdir): os.makedirs(splitdir)

		run_split2Chr(bamfile, splitdir)

	print '\n********All Pre-process Finished********\n'
