
library(mmap)
library(hash)

parseList=function(fn, n) {
  h=new.env(hash=TRUE, parent=emptyenv(), size=n)
  fp=file(fn, open="r")
  ls=readLines(fp)
  for (i in 1:length(ls)) {
    h[[ls[[i]]]]=i
  }
  return (h)
}

# this is a poor man's object.  It returns a named list of methods,
# and internally holds the state those methods need.

getfuncs=function(datfile, rowfile, colfile) {

  fp=file(datfile, "rb")
  seek(fp, -8, "end")
  nr=readBin(fp, integer())
  nc=readBin(fp, integer())

  mp=mmap(file=datfile, mode=real32())

  rows=parseList(rowfile, 2000000)	
  cols=parseList(colfile, 2000000)

  if (nr != length(rows) || nc != length(cols)) {
    cat("Error!! Something is out of whack. I expected the binary file to have", length(rows), "rows and ", length(cols),
	"cols, but it had ", nr, "rows and ", nc, "cols.")
	quit()
  }
  reclen=length(cols)

  getVal = function(row, col) {
    ridx=rows[[row]]
    cidx=cols[[col]]

    idx=1+(ridx-1)*reclen+(cidx-1) # isn't 1-based indexing fun?
    mp[idx]
  }

  getPos = function(row, col) {
    ridx=rows[[row]]
    cidx=cols[[col]]

    idx=1+(ridx-1)*reclen+(cidx-1) # isn't 1-based indexing fun?
    idx
  }

  getRow = function(row) {
    idx=rows[[row]]
    s=1+(as.double(idx)-1)*as.double(reclen)
    e=as.double(idx)*as.double(reclen)
    print (c("e", idx, reclen, e))
    mp[s:e]
  }

  list(getVal=getVal, getPos=getPos, getRow=getRow)
}

# some examples of how to use this
# create api's for both representations
api_snp_gene=getfuncs('combined.dat', 'snps.txt', 'genes.txt')
api_gene_snp=getfuncs('transpose.dat', 'genes.txt', 'snps.txt')

# call a method
api_snp_gene$getVal('rs3094315', 'ILMN_12601_2480184')
api_gene_snp$getVal('ILMN_12601_2480184', 'rs3094315')

length(api_snp_gene$getRow('rs3094315'))
length(api_gene_snp$getRow('ILMN_12601_2480184'))

snps=readLines('snps.txt')
genes=readLines('genes.txt')

for (i in 1:20) {
  snp=sample(snps, 1) 
  gene=sample(genes, 1)
  t1=api_snp_gene$getVal(snp, gene)
  t2=api_gene_snp$getVal(gene, snp)
  if (t1!=t2) {
    print (c(snp, gene, t1, t2))
  }
}

