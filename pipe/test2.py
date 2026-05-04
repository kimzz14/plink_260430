import math

fin = open('gwas/gwas.HD.STD.glm.linear')
fout = open('gwas/filtered.list', 'w')

legend_LIST = fin.readline().rstrip('\n').split('\t')
fout.write('\t'.join(legend_LIST) + '\n')

for line in fin:
    data_LIST = line.rstrip('\n').split('\t')

    if data_LIST[14] == 'NA':continue
    pValue = float(data_LIST[14])
    
    mLogP = -math.log10(pValue)
    
    if mLogP < 4:
        continue

    fout.write(line)

fout.close()