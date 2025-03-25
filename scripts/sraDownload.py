import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf

sralist, outputdir = [sys.argv[1],sys.argv[2]]
if not os.path.exists(outputdir):os.makedirs(outputdir)
# samples correspond to Het_1, Het_2, Imm_1, Imm_2
#with open(sralist, 'r') as f:
#	sra_numbers = f.readlines().split('\n')
#sra_numbers = ["SRR2121685", "SRR2121686", "SRR2121687", "SRR2121688"]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
prefetch = 'prefetch --option-file ' + sralist
print 'Currently downloading: ' + sralist
cf.run_cmd(prefetch)
#for sra_id in sra_numbers:
#	print ("Currently downloading: " + sra_id)
#	prefetch = "prefetch " + sra_id + ' -O ' + outputdir
#	cf.run_cmd(prefetch)
#    print ("The command used was: " + prefetch)
#    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
	print ("Generating fastq for: " + sra_id)
	fastq_dump = "fastq-dump --outdir fastq  --split-files " + outputdir +" sra_id " + ".sra"
#	print ("The command used was: " + fastq_dump)
#    subprocess.call(fastq_dump, shell=True)
	cf.run_cmd(fastq_dump)
