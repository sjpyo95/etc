##source activate py2
import sys, os
sys.path.insert(0, "/home/sunyme95/scripts/modules/")
import argparse
import config as cf
import time

threadN = 14
jobN = 20
genomefile = '/mnt/mone/PMI/CH/Reference/Genome/genome.fa'
#ponfile = '/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/PON/Mutect2-exome-panel.vcf'
ponfile = '/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/PON/MG_KOR.vcf'
germfile = '/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz'
bialfile = '/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/gnomad/biallelic.gnomad.vcf'
#germfile = '/mnt/mone/PMI/CH/02.Variant_Calling/sunyme95/gnomad/af-only-gnomad.raw.sites.vcf'
bedfile = '/mnt/mone/PMI/CH/Reference/On_Target_Bed/Targets-XGEN.66F814ECD1264A978BE2D43A40BDAB94.bed'
imgfile = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/genome.fa.img'
#queue = cf.queue(7)

###Run Mutect2 Calling Variants###
def runMutect(infile, outfile, f1r2file, sample, queue):
	if bedfile == '':
		mutect = 'gatk Mutect2 -R ' + genomefile + ' -I ' + infile +' --germline-resource ' + germfile + ' --panel-ofnormals ' + ponfile + ' -O ' + outfile + ' -RF OverclippedReadFilter --f1r2-tar-gz ' + f1r2file
	else:
#		mutect = 'gatk Mutect2 -R ' + genomefile + ' -I ' + infile + ' --germline-resource ' + germfile + ' --panel-of-normals ' + ponfile + '  -O ' + outfile + ' -L ' + bedfile + ' -RF OverclippedReadFilter --genotype-germline-sites true --genotype-pon-sites true --f1r2-tar-gz ' + f1r2file
		mutect = 'gatk Mutect2 -R ' + genomefile + ' -I ' + infile + ' --germline-resource ' + germfile + ' --panel-of-normals ' + ponfile + '  -O ' + outfile + ' -L ' + bedfile + ' --f1r2-tar-gz ' + f1r2file
#	print '\n'+mutect
#	cf.run_cmd(mutect)
	cf.qsub_execute('Mutect2.'+sample, queue, mutect, threadN, '/'.join(outfile.split('/')[:-1]))
#	cf.qsub_time(jobN, 2)

###Learn the proior probability of read orientation artifact###
def learnReadOrientationModel(f1r2file, sample, queue):
	outfile = f1r2file.split('.f1r2.')[0] + '.read-orientation-model.tar.gz'
	lrom = 'gatk LearnReadOrientationModel -I ' + f1r2file + ' -O ' + outfile
#	print lrom
	cf.qsub_execute('LROM.'+sample, queue, lrom, threadN, '/'.join(f1r2file.split('/')[:-1]))
#	cf.qsub_time(jobN, 2, 'LROM')
#	cf.run_cmd(lrom)
	return outfile

###Filter somatic SNVs and indels###
def filterMutectCalls(vcffile, romfile, outfile, sample,queue):
	filtMC = 'gatk FilterMutectCalls -R ' + genomefile + ' -V ' + vcffile + ' -ob-priors ' + romfile + ' -O ' + outfile
#	print filtMC
	cf.qsub_execute('FMC.'+sample, queue, filtMC, threadN, '/'.join(outfile.split('/')[:-1]))
#	cf.run_cmd(filtMC)
#	cf.qsub_time(jobN, 2, 'FMC')

def FiltAlignArt(bamfile, vcffile, outfile, sample, queue):
	filt_align_art = 'gatk FilterAlignmentArtifacts -R ' + genomefile + ' -V ' + vcffile + ' --bwa-mem-index-image ' + imgfile + ' -I ' + bamfile + ' -O ' + outfile
#	cf.qsub_execute('FAA.'+sample, queue, filt_align_art, threadN, '/'.join(outfile.split('/')[:-1]))

####Left-align indels and trims common bases from indels###
#def leftAlignTrim(infile, outfile, sample, queue):
#	latv = 'gatk LeftAlignAndTrimVariants -R ' + genomefile + ' -V ' + infile + ' -O ' + outfile
##	print latv
	cf.qsub_execute('LATV.'+sample, queue, latv, threadN, '/'.join(outfile.split('/')[:-1]))
##	cf.qsub_time(jobN, 2, 'LATV')
##	cf.run_cmd(latv)
#	return outfile

###Run Mutect2 in multi-samples###
def multiSample(inputdir, outputdir):
	soutputdir = outputdir + '/1.raw/'
	soutputdir2 = outputdir +'/2.f1r2/'
	if not os.path.exists(soutputdir):os.makedirs(soutputdir)
	if not os.path.exists(soutputdir2):os.makedirs(soutputdir2)

	inputfiles = filter(lambda x: '.bam' in x and not '.bai' in x, os.listdir(inputdir))
	f1r2_files = []
	raw_files = []
	for i in range(len(inputfiles)):
		if i % 2 == 1: queue = cf.queue(6)
		else: queue = cf.queue(7)

		infile = inputdir + '/' + inputfiles[i]
		sample = inputfiles[i].split('IDT_')[0]+inputfiles[i].split('_')[1]
		outfile = soutputdir + sample + '.m2.raw.vcf'
		raw_files.append(outfile)
		f1r2file = soutputdir2 + sample + '.f1r2.tar.gz'
		f1r2_files.append(f1r2file)
		runMutect(infile, outfile, f1r2file, sample,queue)

	cf.qsub_time(10, 1, 'Mutect2')

	ori_files = []
	for i in range(len(f1r2_files)):
		if i% 2 == 1: queue = cf.queue(6)
		else: queue = cf.queue(7)
		f1r2_file = f1r2_files[i]
		sample = f1r2_file.split(soutputdir2)[-1].split('.f1r2.tar.gz')[0]
		ori_file = learnReadOrientationModel(f1r2_file, sample, queue)
		ori_files.append(ori_file)

	cf.qsub_time(10,1, 'LROM')

	for i in range(len(raw_files)):
		if i % 2 == 1: queue = cf.queue(6)
		else: queue = cf.queue(7)
		rawfile = raw_files[i]
		sample = rawfile.split(soutputdir)[-1].split('.m2.raw.vcf')[0]
		ori_file = filter(lambda x: sample in x, ori_files)[0]
		
		soutputdir3 = outputdir + '/3.filt/'+sample+'/'
		if not os.path.exists(soutputdir3):os.makedirs(soutputdir3)

		filt_file = soutputdir3+ sample + '.m2.f.vcf'
		filterMutectCalls(rawfile, ori_file, filt_file, sample, queue)
	
	cf.qsub_time(10,1, 'FMC')

	for i in range(len(inputfiles)):
		if i % 2 == 1: queue = cf.queue(6)
		else: queue = cf.queue(7)

		bamfile = inputdir + '/' + inputfiles[i]
		sample = inputfiles[i].split('IDT_')[0]+inputfiles[i].split('_')[1]
		soutputdir4 = outputdir + '/4.algin_filt/'+sample+'/'
		if not os.path.exists(soutputdir4):os.makedirs(soutputdir4)
		vcffile = soutputdir3 + filter(lambda x: '.vcf' in x and not '.idx' in x and not '.tbi' in x, os.listdir(soutputdir3))[0]
		outfile = soutputdir4 + sample + '.m2.f.f.vcf'
		FiltAlignArt(bamfile, vcffile, outfile, sample, queue)

#		pileupt = soutputdir3 + sample + '.pileups.table'
#		pileupts.append(pileupt)
#		getPileup(infile, pileupt, sample, queue)
#	cf.qsub_time(1,1,'getPileup')
#	soutputdir4 = outputdir + '/2_2.contam_table/'
#	if not os.path.exists(soutputdir4): os.makedirs(soutputdir4)
#	for i in range(len(pileupts)):
#		pileupt = pileupts[i]
#		sample = pileupt.split(soutputdir3)[-1].split('.pileups.')[0]
#		contamfile = soutputdir4 + sample + '.contamination.table'
#		calContam(pileupt, contamfile, sample, queue)

###Run Mutect2 in single-sample###
def singleSample(infile, outputdir):
	soutputdir = outputdir + '/1.raw/'
	soutputdir2 = outputdir +'/2.f1r2/'
	if not os.path.exists(soutputdir):os.makedirs(soutputdir)
	if not os.path.exists(soutputdir2):os.makedirs(soutputdir2)
	sample = infile.split('/')[-1].split('IDT_')[0]+infile.split('/')[-1].split('_')[1]
	outfile = soutputdir + '/' + sample + '.m2.raw.vcf'
	f1r2file = soutputdir2 + sample + '.f1r2.tar.gz'
	runMutect(infile, outfile, f1r2file, sample)

###Filtration steps###
def afMutect(infile, outputdir):
	queue = cf.queue(6)
	soutputdir = outputdir + '/1.raw/'
	soutputdir2 = outputdir +'/2.f1r2/'
	sample = infile.split('/')[-1].split('IDT_')[0]+infile.split('/')[-1].split('_')[1]
#	sample = infile.split('/')[-1].split('.m2')[0]
	outfile = soutputdir + '/' + sample + '.m2.raw.vcf'
	f1r2file = soutputdir2 + sample + '.f1r2.tar.gz'
	
	contamfile = outputdir + '/2_2.contam_table/' + sample + '.contamination.table'
	
	cf.qsub_time(1,1,'Mutect2')
	
	romfile = learnReadOrientationModel(f1r2file, sample, queue)
	cf.qsub_time(1,1, 'LROM')
	
	soutputdir = outputdir +'/3.filtflag/'
	if not os.path.exists(soutputdir):os.makedirs(soutputdir)
	rawfile = outfile
	outfile = soutputdir + sample + '.m2.f.vcf'
	
	filterMutectCalls(rawfile, romfile, contamfile, outfile, sample, queue)
	cf.qsub_time(1,1, 'FMC')
	
#	soutputdir = outputdir + '/4.filtout/'+sample+'/'
#	if not os.path.exists(soutputdir):os.makedirs(soutputdir)
#	infile = outfile
#	outfile = soutputdir + sample + '.m2.f.la.vcf'
	
#	vcf = leftAlignTrim(infile, outfile, sample, queue)
#	cf.qsub_time(1,1, 'LATV')
	
	return vcf


if __name__ == '__main__':
	p = argparse.ArgumentParser(description = "Run Mutect2 for single sample or multiple samples")
	p.add_argument('-s', '--single', help = 'If you wish to run in single sample, give this option', action = 'store_true', dest = 'single')
	p.add_argument('-m', '--multi', help = 'Run Mutect2 in multiple samples in single VCF file (no gatk MergeVcfs)', action = 'store_true', dest = 'multi')
	p.add_argument('-I', '--input', help = 'input directory or file', dest = 'input')
	p.add_argument('-O', '--output', help = 'output directory', dest = 'output')
	arg = p.parse_args()
	input = arg.input; output = arg.output; single = arg.single
	multi = arg.multi
	if not os.path.exists(output): os.makedirs(output)
	if single: 
		singleSample(input,output)
		afMutect(input, output)
	else:
		if multi:
			infile = multiInOne(input,output)
			afMutect(infile,output)
		else:	
			multiSample(input,output)

				
	cf.qsub_time(1,1,'psj')


