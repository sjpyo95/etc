import sys, os
sys.path.insert(0, "/home/sunyme95/scripts/modules/")
import config as cf

def runMosdepth(window, fasta, prefix, bamfile,  queue, thread, fastmode, mapq):
	mosdepth = 'mosdepth --by ' + window + ' -f ' +fasta + ' -t ' + thread + ' --mapq ' + mapq + ' '
	if fastmode:
		mosdepth += ' --fast-mode '
	
	mosdepth += prefix + ' ' + bamfile
	cf.qsub_execute('mosdepth', queue, mosdepth, thread, prefix)

inputdir = sys.argv[1]
outputdir = sys.argv[2]
fasta = sys.argv[3]
window = sys.argv[4]
thread = sys.argv[5]
fastmode = sys.argv[6]
mapq = sys.argv[7]
q = int(sys.argv[8])

if fastmode == 'fastmode': fastmode = True
else: fastmode = False

if not os.path.exists(outputdir): os.makedirs(outputdir)

samples = filter(lambda x: '-' in x, os.listdir(inputdir))


for i in range(len(samples)):
	sample = samples[i]
	sinputdir = inputdir + sample + '/'
	bamfile = sinputdir + filter(lambda x: '.bam' in x and not '.bai' in x or '.cram' in x and  not '.crai' in x, os.listdir(sinputdir))[0]
	soutputdir = outputdir + sample + '/'
	if not os.path.exists(soutputdir): os.makedirs(soutputdir)
	prefix = soutputdir + sample
	queue = cf.queue(q)
	runMosdepth(window, fasta, prefix, bamfile, queue, thread, fastmode, mapq)
	

