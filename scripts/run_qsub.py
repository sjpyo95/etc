import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf
import time
import commands

job = sys.argv[1]
q = int(sys.argv[2])
threadN = sys.argv[3]
logdir = sys.argv[4]
cmd = sys.argv[5:]
cmd = ' '.join(cmd)

queue = cf.queue(q)
cf.qsub_execute(job, queue, cmd, threadN, logdir)
