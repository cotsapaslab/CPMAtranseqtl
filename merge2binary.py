import sys, struct, numpy, argparse

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False 

def multi_delete(Lista, idx2remove):
    indexes = sorted(idx2remove, reverse=True)
    for index in indexes:
        del Lista[index]
    return Lista

parser = argparse.ArgumentParser(description='Merges floats in plink assoc files to a binary file')

prefix=sys.argv[1]
suffix=sys.argv[2]
numFiles=sys.argv[3]
binaryFileName=sys.argv[4]
numVarPerFile=sys.argv[5]
kolumn=int(sys.argv[6])
if len(sys.argv)==8:
	removeFile=sys.argv[7]
else:
	removeFile=False

fps=[prefix+'%d'%i+suffix for i in range(1,(int(numFiles)+1))]

removeIdx=[]
if removeFile!=False:
	for line in open(removeFile,'r'):
		removeIdx.append(int(line.strip('\n'))-1) # because these an 1-based indecies, and python is 0-based. 


ofp = open(binaryFileName, 'wb')

for fp in fps:
	vals=[]
	fp=open(fp,'r')
	lines = fp.readlines()
	if not is_number(lines[0].split()[0]):
		lines=lines[1:]
	for l in lines:
		if not l:
			raise StopIteration
		val = l.rstrip().split()[kolumn-1].strip(' ').strip('\n')
		if val=='NA':
			vals.append(val)
		else:
			vals.append(float(val))
	if len(vals)!=int(numVarPerFile):
		print(len(vals),fp)
		break
	vals=multi_delete(vals,removeIdx)
	buf=struct.pack("%df"%len(vals), *vals)
	ofp.write(buf)

ofp.close()
