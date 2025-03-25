import sys

def parse_mpileup(mpileup_file):
    with open(mpileup_file) as f:
        for line in f:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = line[1]
            ref = line[2]
            alt = []
            depth = []
            for i in range(3, len(line), 3):
                if line[i] != '*':
                    alt.append(line[i])
                    depth.append(line[i+1])
            alt_depth = ';'.join([f'{a}:{d}' for a,d in zip(alt, depth)])
            ref_depth = f'REF({ref}):{line[3]};'
            yield f'{chrom}\t{pos}\t{ref}\t{",".join(alt)}\t{ref_depth}{alt_depth}'

if __name__ == '__main__':
    mpileup_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(output_file, 'w') as f:
        f.write('CHROM\tPOS\tREF\tALT\tDEPTH\n')
        for row in parse_mpileup(mpileup_file):
            f.write(row + '\n')

