
import struct, mmap, random, sys

class API:
    def __init__(self, datfile, rowfile, colfile):
        self.datfile=datfile
        self.row_list, self.row_dict=self.parseList(rowfile)
        self.col_list, self.col_dict=self.parseList(colfile)
        self.datfp=open(datfile, 'rb')
        self.map=mmap.mmap(self.datfp.fileno(), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)

        chk_rows, chk_cols=struct.unpack('ii', self.map[-8:])
        self.numcols=len(self.col_dict)
        self.numrows=len(self.row_dict)
        if (chk_rows!=self.numrows or chk_cols!=self.numcols):
            print '''Error!! Something is out of whack. I expected the binary file to have %d rows and %d cols,
            but it had %d rows and %d cols!''' % (self.numrows, self.numcols, chk_rows, chk_cols)
            sys.exit(-1)
        if self.map.size() != self.numcols*self.numrows*4+8:
            print 'Error!! binary file has len %d, I expected %d' % (self.map.size(), self.numcols*self.numrows*4+8)
            sys.exit(-1)
            
        self.rec_len=self.numcols*4
        
    def getRow(self, row):
        idx=self.row_dict[row]
        return struct.unpack('%df'%self.numcols, self.map[self.rec_len*idx:self.rec_len*(idx+1)])

    def getRowAsDict(self, row):
        vals=self.getRow(row)
        d={}
        for k, idx in self.col_dict.iteritems():
            d[k]=vals[idx]
        return d
    
    def getVal(self, row, col):
        idx=self.row_dict[row]*self.rec_len+self.col_dict[col]*4
        return struct.unpack('f', self.map[idx:idx+4])

    def getValRC(self, ridx, cidx):
        assert(isinstance(ridx, int))
        assert(isinstance(cidx, int))
        idx=ridx*self.rec_len+cidx*4
        return struct.unpack('f', self.map[idx:idx+4])

    def getPos(self, row, col):
        return self.row_dict[row]*self.rec_len+self.col_dict[col]*4

    # utility stuff
    def parseList(self, fn):
        d={}
        l=[]
        fp=open(fn)
        for i, name in enumerate(fp):
            name=name.rstrip()
            d[name]=i
            l.append(name)
        return l, d


if __name__=='__main__':
    api_snp_gene=API('combined.dat', 'snps.txt', 'genes.txt')
    api_gene_snp=API('transpose.dat', 'genes.txt', 'snps.txt')

    # test getVal
    for i in range(20):
        snp=api_snp_gene.row_list[random.randint(0, api_snp_gene.numrows-1)]
        gene=api_snp_gene.col_list[random.randint(0, api_snp_gene.numcols-1)]
        print 'testing %s %s' % (snp, gene)
        v1=api_snp_gene.getVal(snp, gene)
        v2=api_gene_snp.getVal(gene, snp)
        assert(v1==v2)

    # test getRow
    print 'genes %d' % len(api_snp_gene.getRow('rs3094315'))
    print 'snps %d' % len(api_gene_snp.getRow('ILMN_11285_5900056'))

    
