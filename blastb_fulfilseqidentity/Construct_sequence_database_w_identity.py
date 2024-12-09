import numpy as np
import sys
from tqdm import tqdm

for i in range(19):
    if i == 0:
       data = np.loadtxt('blastp%d.txt'%i,dtype=str,usecols=(0,1,2))
    else:
       idata = np.loadtxt('blastp%d.txt'%i,dtype=str,usecols=(0,1,2))
       data = np.vstack((data,idata))
indexExcludeSelf = np.argwhere(data[:, 0]!=data[:, 1])
data = data[indexExcludeSelf[:,0],:] 
del indexExcludeSelf
del idata
sequenceIdentity = data[:, 2].astype(np.float64)
identitySelectArray = np.zeros(len(sequenceIdentity),dtype=int)
IndexOfSelection = np.argwhere(sequenceIdentity>50)
identitySelectArray[IndexOfSelection] = 1
data[:,2] = identitySelectArray.astype(str)

del identitySelectArray

delseqname = []
while '1' in data[:,2]:
    data = data[data[:,0].argsort()]
    uniqueSeq, indices = np.unique(data[:,0],return_index=True)
    splitSeq = np.split(data, indices[1:])
    print('uniqueSeq:', len(uniqueSeq))
    print(len(splitSeq),type(splitSeq))
    sumIdentities = [[int(i),np.sum(val[:,2].astype(np.int64))] for i,val in enumerate(splitSeq) if '1' in val[:,2]]
    
    sumIdentities = np.array(sumIdentities,dtype=int)
    minIdentities = np.min(sumIdentities[:,1])
    print('minIdentities:',minIdentities)
    indexMins = np.argwhere(sumIdentities[:,1] == minIdentities) 
    aimChains = [splitSeq[int(sumIdentities[idx,0])] for idx in indexMins]
    
    nametoremove = []
    for aimchain in tqdm(aimChains):
        if aimchain[0,0] in delseqname:
                continue
        indexNotFulfil = np.argwhere(aimchain[:, 2] == '1')
        nameSeqNotFulfil = aimchain[indexNotFulfil[:,0], 1]
        delseqname.extend(nameSeqNotFulfil)

        nametoremove.extend(nameSeqNotFulfil)
    
    isindex = np.isin(data, nametoremove)
    row_to_keep = np.all(~isindex, axis=1)
    data = np.compress(row_to_keep, data, axis=0)
np.savetxt('allsequenceIdentity50.dat',data,fmt='%s')
np.save('alldelsequence.npy', delseqname)





