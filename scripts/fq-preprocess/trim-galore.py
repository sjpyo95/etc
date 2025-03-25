import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf
import time
import commands
import multiprocessing as mp

def runTrimgalore(fastqPairs, outputdir, adapt, trimsize, tN, clip):
#	fastqPairs = datapair[0]; outputdir = datapair[1]; adapt = datapair[2]; trimsize = datapair[3]; tN = datapair[4]; clip = datapair[5]
	fastq1 = filter(lambda x: '_R1_' in x, fastqPairs)[0]
	fastq2 = filter(lambda x: '_R2_' in x, fastqPairs)[0]
	cmd = 'trim_galore ' + adapt + ' --paired -j ' + tN + ' ' + clip + ' -o ' + outputdir + ' ' + fastq1 + ' '+fastq2
#	commands.getoutput(cmd)
	
	return cmd
	

def main(inputdir, outputdir, adapt, trim5size, threadN, jobN, queueN):
	if adapt == 'illumina' : adapt = '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
	if int(trim5size) == 0 : clip = ' '
	else : clip = ' --clip_R1 ' + str(trim5size) + ' --clip_R2 ' + str(trim5size)

	fastqfiles = [inputdir + '/' + x for x in filter(lambda x: '.fastq' in x or '.fq' in x and not '.txt' in x and not '_val_' in x, os.listdir(inputdir))]
	samplesNames = [x.split('/')[-1].split('_R')[0] for x in fastqfiles]
	samples = list(set(samplesNames))
	
	queue = cf.queue(queueN)
	for i in range(len(samples)):
		if i+1%5==0: print i+1,'/', len(samples)
		sample = samples[i]
		fastqs = filter(lambda x: sample in x, fastqfiles)
#		soutputdir = outputdir+ '/' + sample
#		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
		tgCMD = runTrimgalore(fastqs, outputdir, adapt, trim5size, threadN, clip)
		cf.qsub_execute('TG-'+sample, queue, tgCMD, threadN, outputdir)
		cf.qsub_time(jobN, 10)



if __name__ == '__main__':
	inputdir = sys.argv[1]
	outputdir = sys.argv[2]
	if not os.path.exists(outputdir): os.makedirs(outputdir)
	adapt = sys.argv[3]
	trim5size = sys.argv[4] #e.g. 0, 12
	threadN = sys.argv[5]
	jobN = int(sys.argv[6])
	queueN = int(sys.argv[7])

	main(inputdir, outputdir, adapt, trim5size, threadN, jobN, queueN)


