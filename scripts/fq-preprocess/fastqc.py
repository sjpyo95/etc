import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf
import time
from itertools import repeat
import subprocess
#cmd1 = 'export PATH=/home/sunyme95/anaconda3/bin:$PATH'
#cmd2 = 'source activate CH'
#os.system(cmd1)
#os.system(cmd2)

def runFastqc(fastq, outputdir, threadN):
#def runFastqc(datapair):
#	print fastq
#	fastq = datapair[0]; outputdir = datapair[1]; tN = datapair[2]
	fastqcRun ='fastqc' + ' -t ' + threadN + ' -o ' + outputdir + ' ' + fastq
	
#	print fastqcRun
	subprocess.run(fastqcRun, shell=True, check=True)
	return fastqcRun

def main(inputdir, outputdir, threadN):#, jobN, queueN):
	fastqs = [inputdir + '/' + x for x in  filter(lambda x: '.txt' not in x and 'log' not in x and '.fastq' in x or '.fq' in x, os.listdir(inputdir))]
#	data_pairs = zip(fastqs, repeat(outputdir), repeat(threadN))
#	for i in range(0,len(data_pairs), 4):
#	dps = data_pairs[i:i+4]
#	pool = mp.Pool(processes= int(multiN))
#	for i in range(len(data_pairs)):
#		dp = data_pairs[i]
#		sample = dp[0].split('/')[-1].split('.fastq')[0]
#		queue = cf.queue(queueN)
		
#		cf.qsub_time(jobN,2)
#		cf.qsub_execute('fastqc.'+sample, queue, runFastqc(dp), str(int(threadN)+2), outputdir)
	for i in range(len(fastqs)):
		fastq = fastqs[i]
		runFastqc(fastq, outputdir, threadN)
#		if i == 2: exit()
if __name__ == '__main__':

	inputdir = sys.argv[1]
	outputdir = sys.argv[2]
	threadN = sys.argv[3]
#	jobN = int(sys.argv[4])
#	queueN = int(sys.argv[5])
	main(inputdir, outputdir, threadN)#, jobN, queueN)
