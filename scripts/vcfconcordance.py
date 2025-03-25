import vcf
import os, sys
import time
def time_cal(start_time, end_time):
	runtime = end_time - start_time
	hrs = int(runtime // 3600)
	mins = int((runtime % 3600) // 60)
	secs = int(runtime % 60)
	output_string = f"{hrs}hr {mins}min {secs}sec"
	return output_string

def compare_positions_count(vcf1_pos, vcf2_pos, outputdir):
	common_positions = set(vcf1_pos) & set(vcf2_pos)
	vcf1_only = set(vcf1_pos) - common_positions
	vcf2_only = set(vcf2_pos) - common_positions
	with open(outputdir + '/Concordance.txt', 'w') as f:
		f.write('Common Positions\t'+str(len(common_positions))+'\n')
		f.write('Illumina\t'+str(len(vcf1_only))+'\n')
		f.write('Ultima\t'+str(len(vcf2_only)))


start_time = time.time()

print("Start Parsing VCF files...")
vcf1_file = sys.argv[1]
vcf2_file = sys.argv[2]
outputdir = sys.argv[3]
if not os.path.exists(outputdir): os.makedirs(outputdir)

vcf1 = vcf.Reader(open(vcf1_file, 'rb'))
vcf2 = vcf.Reader(open(vcf2_file, 'rb'))

vcf1_pos = [(x.CHROM, x.POS) for x in vcf1]
vcf2_pos = [(x.CHROM, x.POS) for x in vcf2]
print(vcf1_file+':',len(vcf1_pos), "Positions")
print(vcf2_file+':', len(vcf2_pos), "Positions")
print(time_cal(start_time, time.time()), "\nDone Parsing\n")

print("Calculating Concordance...")

compare_positions_count(vcf1_pos, vcf2_pos, outputdir)
print(time_cal(start_time, time.time()), "\nDone")
