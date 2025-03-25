import sys, os
import vcf
#from collections import defaultdict
import time
def time_cal(start_time, end_time):
	runtime = end_time - start_time
	hrs = int(runtime // 3600)
	mins = int((runtime % 3600) // 60)
	secs = int(runtime % 60)
	output_string = f"{hrs}hr {mins}min {secs}sec"
	return output_string

def write_table(filename, data):
	with open(filename, 'w') as f:
		f.write("CHROM\tPOS\tREF\tALT1\tALT2\tDP1\tDP2\tDifference\n")
		for row in data:
			f.write('\t'.join(map(str, row))+'\n')

def main(vcf1_file, vcf2_file, outputdir):
	start_time = time.time()
	print('Reading VCF Files...')
	vcf1 = vcf.Reader(open(vcf1_file, 'r'))
	vcf2 = vcf.Reader(open(vcf2_file, 'r'))
	
	readvcftime = time.time()
	outtime = time_cal(start_time, readvcftime)
	
	print('Done\n')
	common_nonrefs = []
	comm_vcf1_alt = []
	comm_vcf2_alt = []
	ref_unmatch = []
	alt_unmatch = []
	match = []
	vcf1_only = []
	vcf2_only = []
	vcf1count = 0
	vcf2count = 0
	print("Calculating Each positions...")

	for record in vcf2:
		vcf2count += 1
		vcf2_records[(record.CHROM, record.POS)] = record
	vcf2_only = vcf2_records.keys()

	print(vcf2_file, len(vcf2_records.keys()), "Positions...")
	parsetime = time.time()
	outtime = time_cal(readvcftime, parsetime)
	print(outtime,'...\n')

	for record in vcf1:
		vcf1count += 1
		key = (record.CHROM, record.POS)
		dp1 = sum(record.samples[0]['AD'])
		dp2 = 'NA'
		if key in vcf2_records.keys():
			matching_record = vcf2_records[key]
			dp2 = sum(record.samples[0]['AD'])

			if record.REF == record2.REF:
				common_all.append((record.CHROM, record.POS, record.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1 - dp2))
#				if None in record.ALT and len(record2.ALT) == 1:
#					common_nonrefs.append((record.CHROM, record.POS, record.REF, '.', str(record2.ALT), dp1, dp2, dp1 - dp2))
#				elif None in record.ALT:
#					comm_vcf2_alt.append((record.CHROM, record.POS, record.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1 - dp2))
#				elif len(record2.ALT) == 1:
#					comm_vcf1_alt.append((record.CHROM, record.POS, record.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1-dp2))
#				if record.ALT != record2.ALT:
#					alt_unmatch.append((record.CHROM, record.POS, record.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1 - dp2))
#				if record.ALT == record2.ALT:
#					match.append((record.CHROM, record.POS, record.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1 - dp2))

		else:
			ref_unmatch.append((record.CHROM, record.POS, record.REF, record2.REF, str(record.ALT), str(record2.ALT), dp1, dp2, dp1 - dp2))
				found_match = True
				break

		if found_match:
			continue

		if not matching_records:
			vcf1_only.append((record.CHROM, record.POS, record.REF, str(record.ALT), 'NA', dp1, 'NA', dp1))

	common_pos = [(x[0], x[1]) for x in common_nonrefs+comm_vcf1_alt+comm_vcf2_alt]

	if key in vcf2_records.keys():
		if key in [(x[0], x[1]) for x in common
		vcf2_only.remove(key)
	
	caltime = time.time()
	outtime = time_cal(parsetime, caltime)
	print('Done Calculating..\n', outtime,'\n')

	with open(outputdir+'stat.txt', 'w') as f:

		f.write("Total Illumina\t"+str(vcf1count)+'\n')
		f.write("Total Ultima\t"+str(vcf2count)+'\n')
	#print("common position\t"+str(len(common_all)))
		f.write("common reference\t"+str(len(common_nonrefs))+'\n')
		f.write("Matched Alternative positions\t"+str(len(match))+'\n')
		f.write("Unmatched Alternative positions\t"+str(len(alt_unmatch))+'\n')
		f.write("Variant only in Illumina\t"+str(len(comm_vcf1_alt))+'\n')
		f.write("Variant only in Ultima\t"+str(len(comm_vcf2_alt))+'\n')
		f.write("Unmatched reference base\t"+str(len(ref_unmatch))+'\n')
		f.write("Illumina only\t"+str(len(vcf1_only))+'\n')
		f.write("Ultima only\t"+str(len(vcf2_only)))

	write_table(outputdir+"/common_nonrefs.txt", common_nonrefs)
	write_table(outputdir+"/matched_alternative_positions.txt", match)
	write_table(outputdir+"/unmatched_alternative_positions.txt", alt_unmatch)
#	write_table(outputdir+"/common_alts.txt", common_alts)
	write_table(outputdir+"/common_Illumina_variant.txt", comm_vcf1_alt)
	write_table(outputdir+"/common_Ultima_variant.txt", comm_vcf2_alt)
	write_table(outputdir+"/Illumina_only.txt", vcf1_only)
	
	with open(outputdir+"/Ultima_only.txt", 'w') as f:
		f.write("CHROM\tPOS\tREF\tALT1\tALT2\tDP1\tDP2\tDifference\n")
		for key in vcf2_only:
			records = vcf2_records.get(key, [])
			for record in records:
				f.write('\t'.join(map(str, (record.CHROM, record.POS, record.REF, 'NA', str(record.ALT), sum(record.samples[0]['AD']), -sum(record.samples[0]['AD'])+'\n'))))
	
	endtime = time.time()
	outtime = time_cal(start_time, endtime)
	print(f"Total running time:{outtime}")

if __name__ == "__main__":

	vcf1_file = sys.argv[1]
	vcf2_file = sys.argv[2]
	outputdir = sys.argv[3]
	if not os.path.exists(outputdir): os.makedirs(outputdir)

	main(vcf1_file, vcf2_file, outputdir)
