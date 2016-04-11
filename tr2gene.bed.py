import sys

krom=set()
for i in range(1,23):
	krom.add(''.join(['chr',str(i)]))

krom.add('chrX')
krom.add('chrY')

genes={}

for line in open(sys.argv[1],'r'):
	if not line.startswith('#bin'):
		if line.split('\t')[2] in krom:
			if line.split('\t')[12] not in genes.keys():
				genes[line.split('\t')[12]]={'strand':set(),'chr':set(),'start':set(), 'end':set()}
			genes[line.split('\t')[12]]['start'].add(line.split('\t')[4])
			genes[line.split('\t')[12]]['end'].add(line.split('\t')[5])
			genes[line.split('\t')[12]]['strand'].add(line.split('\t')[3])
			genes[line.split('\t')[12]]['chr'].add(line.split('\t')[2])

outfile=open('.'.join([sys.argv[1],'hgnc']),'w')

for gene in genes.keys():
	if len(genes[gene]['chr'])==1:
		outfile.writelines('\t'.join([min(genes[gene]['chr']),str(min([int(x) for x in genes[gene]['start']])), str(max([int(x) for x in genes[gene]['end']])),gene])+'\n')

outfile.close()
