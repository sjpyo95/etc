#!/usr/bin/python3
import sys, os
from subprocess import Popen, PIPE, STDOUT
import subprocess
import time
import re
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial

hm_chrs = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

def convert_time(sec):
	sec = sec % (24*3600)
	hour = sec // 3600
	sec %= 3600
	min = sec // 60
	sec %= 60
	return "%02d:%02d:%02d" % (hour, min, sec)

def sortChrPos(chromosome_position_list):
	return sorted(chromosome_position_list, key=lambda x: (int(x.split(":")[0]) if x.split(":")[0].isdigit() else float('inf'), int(x.split(":")[1])))
####DATA BASE####

#hg19 VCF Files
NARD2 = '/mnt/mone/PMI/CH/01.Alignment/khyojoo01/001-REFGENOME/NARD2_MAF.hg19.sorted.vcf.gz'
dbsnp = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/gatk4_compatible/dbsnp_132_b37.leftAligned.vcf.gz'
indel_1000g = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/gatk4_compatible/1000G_phase1.indels.hg19.sites.fixed.vcf.gz'
Mills_1000g = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/gatk4_compatible/Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf.gz'
hapmap = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/hapmap_3.3.b37.vcf.gz'
gnomad_exom  = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/gnomad.exomes.r2.1.1.sites.vcf.bgz'
omni = '/mnt/mone/PMI/CH/Reference/Genome/bwa-mem/gatk4_compatible/1000G_phase1.indels.hg19.sites.fixed.vcf.gz'

#def splitThreads(totalThread, jobN):
#	totalThread/jobN

def queue(i):
	if i % 6 == 5: queue = 'bd1.q@bdcm02'
	elif i % 6 == 4: queue = 'bd1.q@bdcm02'
	elif i % 6 == 3: queue = 'bd1.q@bdcm04'
	elif i % 6 == 2: queue = 'bd1.q@bdcm04'
	elif i % 6 == 1: queue = 'bd1.q@bdcm05'
#	elif i % 7 == 1: queue = 'bd1.q@bdcm08'
	else: queue = 'bd1.q@bdcm05'
#	if i == 1 : queue = 'bd1.q@bdcm01'
#	if i == 2: queue = 'bd1.q@bdcm02'
#	elif i == 4 : queue = 'bd1.q@bdcm04'
#	elif i == 5: queue = 'bd1.q@bdcm05'
#	elif i == 6 : queue ='bd1.q@bdcm06'
#	elif i == 7 : queue ='bd1.q@bdcm07'
#	elif i == 8 : queue ='bd1.q@bdcm08'
#	elif i == 0 : queue = 'bd1.q'
#	else : print 'Error: Check your queue number!'; exit()
	return queue

def qsub_execute(job, queue, cmd, threadN, logdir):
#	if not os.path.exists(logdir+'/logs/') : os.makedirs(logdir+'/logs/')
	qsub = 'qsub -q ' + queue + ' -pe pePAC ' + str(threadN) + ' -cwd  -V ' + ' -o ' + logdir +  '.out -e ' + logdir + '.err.out -b y -N ' + 'psj.' + job + ' "' + cmd +'"'
	print( '\n' + qsub + '\n')
	p = Popen(qsub, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()

def qsub_time(jobN, t, job):
#	while commands.getoutput('qstat').count(job) >= jobN: time.sleep(t)
	while True:
		qstat_output = subprocess.check_output(["qstat"]).decode()
		job_count = qstat_output.count(job)
		if job_count >= jobN:
			time.sleep(t)
		else:
			break

def run_cmd(cmd):
	print( '\n' + cmd + '\n')
	p = Popen(cmd, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
	p.wait()


def abortable_worker(func, *args, **kwargs):
	timeout = kwargs.get('timeout', None)
	coreN = kwargs.get('threadN',None)
	p = ThreadPool(coreN)
	res = p.apply_async(func, args = args)
	try:
		out = res.get(timeout)
		return out

	except mp.TimeoutError:
		print( args, (': Aborting due to timeout'))
		raise

def multi_work(data_pairs, func, coreN, timeoutN):
	pool = mp.Pool(processes = coreN)
	for dp in data_pairs:
		abortable_func = partial(abortable_worker, func, timeout=timeoutN, threadN=coreN)
		pool.apply_async(abortable_func, args = dp)

	pool.close()
	pool.join()
