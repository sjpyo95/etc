import sys, os
import gzip
sys.path.insert(0, '/home/sunyme95/scripts/modules/')
import config as cf
import commands
import time

def getInfo(info):
	infos = info.split(';')
	ac = int(infos[0]); af = float(infos[1]); an = infos[2]; baseQRankSum = infos[3]
	return ac,af

def changeScore2Nucl(score, ref, alt):
	scores = score.split('/')
	rscore = int(scores[0]); asocre = int(scores[1])
	if ascore == 0:
		return str(ref) + '/' + str(ref)
	elif rscore > 0 :
		return ref + '/' + alt


def parseVCF(gzipfile):
	infile = gzip.open(gzipfile, 'rb')
	lines = infile.readlines(); infile.close()
	samdic = dict()
	for i in range(len(lines)):
		line = lines[i].strip()
		if line.startswith('#CHROM'):
			tmp = line.split('\t')
			samples = tmp[9:]
		elif not line.startswith('#'):
			tmp = line.split('\t')
			chr = tmp[0]; pos = tmp[1]; id = tmp[2]
			refnucl = tmp[3]; altnucl = tmp[4]; qual = tmp[5]; filt = tmp[6]; info  = tmp[7]; form = tmp[8]
			if not altnucl == '<NON_REF>': continue
			altnucl = refnucl
#			if filt != 'PASS': continue

#			ac, af = getInfo(info)
			cols = (chr,pos,id,refnucl,altnucl,qual,filt,info,form)
#			region = chr+':'+pos
			if not samdic.has_key(cols):
				samdic[cols] = dict()
				
			sampleScores = tmp[9:]
			for t in range(len(samples)):
				sample = samples[t]; score = sampleScores[t]
				score_tmp = score.split(':')
				score = '0/0:'+':'.join(score_tmp[1:])

				samdic[cols][sample] = score
	return samdic, samples



def main(vcf, samples, outfilename):
	outgzip = gzip.open(outfilename, 'wb')
	outgzip.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' + '\t'.join(samples)+'\n')
	for cols in vcf.keys():
		scores = vcf[cols]
		outgzip.write('\t'.join(list(cols)))
		for sample in samples:
			score = scores[sample]
			outgzip.write('\t'+score)
		outgzip.write('\n')
	outgzip.close()


if __name__ == '__main__':
	VCFfile = sys.argv[1]
	outfile = sys.argv[2]
	vcf,samples = parseVCF(VCFfile)
	main(vcf,samples,outfile)


