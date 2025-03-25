import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config
import time
import commands


def runBWA(fastq1, fastq2, refGenome, outputdir, tN, sample):
#	fastq1 = dp[0]; fastq2 = dp[1]; refGenome = dp[2]; outputdir = dp[3]; tN = dp[4]; sample = dp[5]
	if not os.path.exists(outputdir): os.makedirs(outputdir)

	bwamemCMD = 'bwa mem -M -t ' + tN + ' ' + refGenome + ' ' + fastq1 + ' ' + fastq2  + ' > ' + outputdir + '/' + sample + '.sam'
	return bwamemCMD


def main(inputdir, outputdir, refGenomefile, tN, jobN, queueN):
	fastqfiles = filter(lambda x: '.txt' not in x and 'log' not in x and '.fq' in x, os.listdir(inputdir))
	samplesNames = [x.split('/')[-1].split('_R')[0] for x in fastqfiles]
	samples = list(set(samplesNames))
	queue = config.queue(queueN)
	for i in range(len(samples)):
		print i+1,'/',len(samples)
		sample = samples[i]
		fastqs = [inputdir + x for x in filter(lambda x: sample in x, fastqfiles)]
		fastq1 = filter(lambda x: '_R1_' in x, fastqs)[0]
		fastq2 = filter(lambda x: '_R2_' in x, fastqs)[0]
		soutputdir = outputdir + '/' + sample
#		if not os.path.exists(soutputdir): os.makedirs(soutputdir)
		bwaCMD = runBWA(fastq1,fastq2, refGenomefile, soutputdir, tN, sample)
		config.qsub_execute('bwa-'+sample, queue, bwaCMD, tN, outputdir)
		config.qsub_time(jobN, 20)

#	for data in data_pairs:
#	pool = mp.Pool(processes = multiN)
#	p = pool.map(runBWA, data_pairs)

#	pool.terminate()

if __name__== '__main__':
	inputdir = sys.argv[1]
	outputdir = sys.argv[2]
	refGenomefile = sys.argv[3]
	tN = sys.argv[4]
	jobN = int(sys.argv[5])
	queue = int(sys.argv[6])

	
	main(inputdir, outputdir, refGenomefile, tN, jobN, queue)

