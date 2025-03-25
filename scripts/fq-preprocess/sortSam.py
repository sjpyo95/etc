import sys, os
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf

def sortSAM(samfile, sample, outputdir, tN):
	soutputdir1 = outputdir + '/bam/'+sample+'/'
	if not os.path.exists(soutputdir1): os.makedirs(soutputdir1)
	soutputdir2 = outputdir + '/dpmrk/'+sample+'/'
	if not os.path.exists(soutputdir2): os.makedirs(soutputdir2)

	sam2bam = 'samtools view -@ ' + tN + ' -bS ' + samfile + ' | samtools sort -@ ' + tN + ' -n -o ' + soutputdir1+sample+'.bam'
	dpmrk = 'samtools collate -@ ' + tN + ' -O -u ' + soutputdir1+sample+'.bam | samtools fixmate -@ ' + tN + ' -m -u - - | samtools sort -@ ' + tN + ' -u - | samtools markdup -@ ' + tN + ' - ' + soutputdir2+sample+'.dpmrk.bam'
	cf.run_cmd(sam2bam)
	cf.run_cmd(dpmrk)
	
	idx = 'samtools index -@ ' + tN + ' ' + soutputdir2+sample+'.dpmrk.bam'
	cf.run_cmd(idx)
	
if __name__ == '__main__':
	samfile = sys.argv[1]
	outputdir = sys.argv[2]
	sample = sys.argv[3]
	threadN = sys.argv[4]
#	multiN = int(sys.argv[4])
	sortSAM(samfile, sample, outputdir, threadN)
	cf.qsub_time(1,1,'sort')
#	main(inputdir, outputdir, threadN, multi)
