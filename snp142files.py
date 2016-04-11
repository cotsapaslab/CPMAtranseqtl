# A script that takes the complete snp142 data set and produces files that plink can use 
# to update annotation.

krom=set()
for i in range(1,23):
	krom.add(''.join(['chr',str(i)]))

krom.add('chrX')
krom.add('chrY')

outfile=open('snp142_chr.txt','w')
for line in open('snp142.txt','r'):
	if line.split('\t')[1] in krom:
		outfile.writelines('\t'.join([line.split('\t')[4],line.split('\t')[1].strip('chr')])+'\n')

outfile.close()

outfile=open('snp142_bd.txt','w')
for line in open('snp142.txt','r'):
	if line.split('\t')[1] in krom:
		outfile.writelines('\t'.join([line.split('\t')[4],line.split('\t')[2]])+'\n')

outfile.close()
