#### This document outlines the code used to retrieve and create annotation files used in the trans-eQTL study. 
#### Author: Boel Brynedal 
#### Files created and provided along with code:
#### geneAnno.txt
#### gencode20.hg38.hgnc.sorted.bed
####

#### Collect HGNC symbols for illumina probes

# Load illumina probes:
commonGenes=read.table('commonGenes',header=F,as.is=T,stringsAsFactors=F)

library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
geneAnno=getBM(attributes=c('illumina_humanwg_6_v2', 'hgnc_symbol'),filters = 'illumina_humanwg_6_v2', values = commonGenes, mart = ensembl)
geneAnno2=cbind(commonGenes,geneAnno[match(commonGenes,geneAnno[,1]),2])
geneAnno=geneAnno2
colnames(geneAnno)=c('ID','HGNC')

# Some of the resulting HGNC symbols did not match GENCODE names. Since the TF analysis relies on GENCODE annotation we decided to stick to GENCODE v 20 gene symbols. 
# Update some HGNC symbols to GENCODE v20 names.
geneAnno[which(geneAnno[,2]=='CFAP36'),2]='CCDC104'
geneAnno[which(geneAnno[,2]=='DRC7'),2]='CCDC135'
geneAnno[which(geneAnno[,2]=='ARRDC1-AS1'),2]='C9orf37'
geneAnno[which(geneAnno[,2]=='C11orf98'),2]='C11orf48'
geneAnno[which(geneAnno[,2]=='MTERF3'),2]='MTERFD1'
geneAnno[which(geneAnno[,2]=='CFAP45'),2]='CCDC19'
geneAnno[which(geneAnno[,2]=='ZNF271P'),2]='ZNF271'
geneAnno[which(geneAnno[,2]=='NAPRT'),2]='NAPRT1'
geneAnno[which(geneAnno[,2]=='URB1-AS1'),2]='C21orf119'
geneAnno[which(geneAnno[,2]=='CFAP20'),2]='C16orf80'
geneAnno[which(geneAnno[,2]=='BCO1'),2]='BCMO1'
geneAnno[which(geneAnno[,2]=='GATB'),2]='PET112'
geneAnno[which(geneAnno[,2]=='SARAF'),2]='TMEM66'
geneAnno[which(geneAnno[,2]=='CFAP57'),2]='WDR65'

# Missing in GENCODE v20, therefore deleted. 
geneAnno[which(geneAnno[,2]=='RPL37P6'),2]=NA
geneAnno[which(geneAnno[,2]=='HSP90AA2P'),2]=NA
geneAnno[which(geneAnno[,2]=='IMMP1LP1'),2]=NA


save(geneAnno,file='geneAnno.RData')
write.table(geneAnno,file='geneAnno.txt',quote=F,col.names=F,row.names=F)


# Collect GENCODE v20 annotation in BED format from USCS tables. This file is also supplied:
# gencode20.hg38.all
# This file contains multiple transcripts (lines) per gene. The HGNC symbol is in column 12. 
# We wanted a file containing start and end for HGNC genes. We collected the smallest start and the largest end for each gene, and created a new BED file with HGNC names. We only kept annotation for the probe sets we are analyzing, and only on chromosome 1-22, X and Y. Remove genes with transcripts that are not overlapping.   

# Python code:

# Create a library of the genes we are analyzing:
mina={}

for line in open('geneAnno.txt','r'):
	mina[(line.split(' ')[1].strip('\n'))]={'strand':set(),'chr':set(),'start':set(), 'end':set()}

# Only keep genes on chr 1-22, X and Y:
krom=set()
for i in range(1,23):
	krom.add(''.join(['chr',str(i)]))

krom.add('chrX')
krom.add('chrY')

# Collect hg38 coordinates for the genes we are analyzing in mina dictionary
# If a transcript does not overlap with other transcripts of the same gene it is put in nm. 

mm=set()

for line in open('gencode20.hg38.all','r'):
	if not line.startswith('#bin'):
		if line.split('\t')[12] in mina:
			if line.split('\t')[2] in krom:
				if len(mina[line.split('\t')[12]]['start'])>=1:
					if (min([int(x) for x in mina[line.split('\t')[12]]['start']])>int(line.split('\t')[5]) or  max([int(x) for x in mina[line.split('\t')[12]]['end']])<int(line.split('\t')[4])):
						mm.add(line.split('\t')[12])
				mina[line.split('\t')[12]]['start'].add(line.split('\t')[4])
				mina[line.split('\t')[12]]['end'].add(line.split('\t')[5])
				mina[line.split('\t')[12]]['strand'].add(line.split('\t')[3])
				mina[line.split('\t')[12]]['chr'].add(line.split('\t')[2])


outfile=open('gencode20.hg38.hgnc.bed','w')

for gene in mina.keys():
	if len(mina[gene]['chr'])==1 and gene not in mm:
		outfile.writelines('\t'.join([min(mina[gene]['chr']),str(min([int(x) for x in mina[gene]['start']])), str(max([int(x) for x in mina[gene]['end']])),gene])+'\n')

outfile.close()

sort -k1,1 -k2,2n gencode20.hg38.hgnc.bed > gencode20.hg38.hgnc.sorted.bed


