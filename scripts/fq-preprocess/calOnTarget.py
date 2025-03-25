import sys, os
import commands
import time

def comma(num):
	'''Add comma to every 3rd digit. Takes int or float and returns string.'''
	if type(num) == int:
		return '{:,}'.format(num)
	elif type(num) == float:
		return '{:,.2f}'.format(num) # Rounds to 2 decimal places
	else:
		print "Need int or float as input to function comma()!"

def basePairs(reportfile):
	origin_line = commands.getoutput('grep "Total basepairs processed" ' + reportfile)
	originbp = origin_line.split(':')[-1].strip()
	trimbp_line = commands.getoutput('grep "Total written" ' + reportfile)
	trimbp = trimbp_line.split(':')[-1].strip()
#	print 'grep "Total basepairs processed" ' + reportfile
#	print 'grep "Total written (filtered)" ' + reportfile
	print 'Origin bp : ', originbp
	print 'Trimmed  bp : ', trimbp
	
	return originbp, trimbp

def duplicateRate(bamfile):
#	print 'samtools flagstat -@ 10 ' + bamfile
	reportlines = commands.getoutput('samtools flagstat -@ 15 ' + bamfile)
	dup = reportlines.split('\n')[3].split('+')[0].strip()
	mapped = reportlines.split('\n')[4].split('+')[0].strip()
	mapEff = mapped + ' (' + reportlines.split('\n')[4].split('(')[-1].split(':')[0].strip()+')'
	dupRate = round(float(dup)/float(mapped)*100, 2)
	print 'Mapped Efficiency : ', mapped, '('+mapEff+')'
	print 'Duplication rate : ', dup, '('+str(dupRate)+'%)'

	return mapEff, dupRate

#	return dupRate

def TargetCov(bedfile, bamfile, tN):
#	print "awk '{sum+=$3-$2} END {print sum}' " + bedfile
#	print "samtools mpileup -Q 20 -q 20 -I " + bedfile + ' ' + bamfile + " | cut -f 1,2,4 | awk '{sum+=$3} END {print sum}'"
	totalLen = float(commands.getoutput("awk '{sum+=$3-$2} END {print sum}' " + bedfile).strip())

	totalCov = commands.getoutput("samtools mpileup -d 10000 -Q 20 -q 20 -l " + bedfile + ' ' + bamfile + " | cut -f 1,2,4 | awk '{sum+=$3} END {print sum}'")
	totalCov = totalCov.split('\n')[-1]

	print 'Total length : ', totalLen
	print 'Total base count : ', totalCov
	
	targetcov = float(totalCov)/totalLen
	
	print 'Target Coverage : ', str(round(targetcov,2))+ 'x'
#	print 'samtools view -@ ' + tN +  ' -Q 20 -q 20 -c ' + bamfile
#	print 'samtools view -@ ' + tN + ' -Q 20 -q 20 -L '+ bedfile + ' -c ' + bamfile
	totalreads = float(commands.getoutput('samtools view -@ ' + tN + ' -q 20 -c ' + bamfile).split('\n')[-1])
	ontargetreads = float(commands.getoutput('samtools view -@ ' + tN + ' -q 20 -L '+ bedfile + ' -c ' + bamfile).split('\n')[-1])
	onTargetP = ontargetreads/totalreads*100
	offTargetP = 100.0-onTargetP
	
	print 'Total reads : ', totalreads
	print 'On-targeted reads : ', ontargetreads, '('+str(round(onTargetP,2))+'%)'

	return round(targetcov, 2), round(onTargetP, 2), round(offTargetP, 2)

#def writeReport(sampleID, originbp, trimbp, mapEff, dupRate, targetcov, onTargetP, offTargetP, filename):
	

def main(bamdir,bedfile,trim_report_dir, outputfile, threadN):
#	start_time = time.time()
	samples = filter(lambda x: 'log' not in x and '.bam' in x and '.bai' not in x, os.listdir(bamdir))
	infodic = dict()
	outfile = open(outputfile, 'w')
	outfile.write('SampleID\tOriginal_Base_Pairs\tBase_Pairs_After_Trimming\tDuplicate_Rate(%)\tMapping_Efficiency(%)\tTarget_Coverage(x)\tOn-target(%)\tOff-target(%)')

	for i in range(len(samples)):
		start_time = time.time()
		bam = samples[i]
		sample = bam.split('.markdup')[0]
		print '\n' + sample + '\n'
		bamfile = bamdir + '/' + bam
		trim_report = trim_report_dir + '/' + sample + '_R1_001.fastq.gz_trimming_report.txt'
		oribp, trimbp = basePairs(trim_report)
		mapEff, dupRate = duplicateRate(bamfile)
		targetcov, onTargetP, offTargetP = TargetCov(bedfile, bamfile, threadN)
		print '\n' + 'Running Time : ',time.time()-start_time,'\n'
		outfile.write('\n' + sample + '\t' + str(oribp) + '\t'+ str(trimbp) + '\t' + str(dupRate) + '\t' + str(mapEff) + '\t' + str(targetcov) + '\t' + str(onTargetP) + '\t' + str(offTargetP))
	outfile.close()

#		writeReport(sample, str(oribp), str(trimbp), str(mapEff), str(dupRate), str(targetcov), str(onTargetP), str(offTargetP), outputfile)

	
if __name__=='__main__':
	bamdir = sys.argv[1]
	bedfile = sys.argv[2]
	trim_report_dir = sys.argv[3]
	outputfile = sys.argv[4]
	threadN = sys.argv[5]
	main(bamdir,bedfile,trim_report_dir, outputfile, threadN)
