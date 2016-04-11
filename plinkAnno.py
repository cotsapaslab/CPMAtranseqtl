# A script that takes a USCS data dump, tab delimited, and create files plink uses to 
# update base pair and chromosomes annotation.
# Only keeps SNPs on chromosomes 1-23, X & Y
# Usage: python input.txt
# Output: input_bp.txt & input_chr.txt

import sys

krom=set()
for i in range(1,23):
	krom.add(''.join(['chr',str(i)]))

krom.add('chrX')
krom.add('chrY')

outfile=open(''.join([sys.argv[1],str('.txt'),'_chr.txt']),'w')
for line in open(sys.argv[1],'r'):
	if line.split('\t')[1] in krom:
		outfile.writelines('\t'.join([line.split('\t')[4],line.split('\t')[1].strip('chr')])+'\n')

outfile.close()

outfile=open(''.join([sys.argv[1],str('.txt'),'_bp.txt']),'w')
for line in open(sys.argv[1],'r'):
	if line.split('\t')[1] in krom:
		outfile.writelines('\t'.join([line.split('\t')[4],line.split('\t')[2]])+'\n')

outfile.close()
