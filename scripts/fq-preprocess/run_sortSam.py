import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf
import time

inputdir = sys.argv[1]
outputdir = sys.argv[2]
threadN = sys.argv[3]
jobN = int(sys.argv[4])
queueN = int(sys.argv[5])

samples = filter(lambda x: '.txt' not in x and 'log' not in x, os.listdir(inputdir))
if not os.path.exists(outputdir): os.makedirs(outputdir)
queue = cf.queue(queueN)
for s in range(len(samples)):
	print s+1,'/',len(samples)
	sample = samples[s]
	sinputdir = inputdir+'/'+sample+'/'
	samfile = sinputdir + filter(lambda x: '.sam' in x, os.listdir(sinputdir))[0]
	cmd = '/usr/bin/python2.7 /home/sunyme95/scripts/fq-preprocess/sortSam.py ' + samfile + ' ' + outputdir + ' ' + threadN
	cf.qsub_execute('sortSAM-'+sample, queue, cmd, str(int(threadN)+1), outputdir)
	cf.qsub_time(jobN, 20)

