#!/mnt/cfs/med18b/med/renqian/env//anaconda3/bin/python
###################################################################
# File Name: neoantigen_report2data.py
# Author: renqian
# mail: renqian@yucebio.com
# Created Time: 2018年05月17日 星期四 10时01分48秒
#=============================================================
import sys,os
infile,output,sample=sys.argv[1:]
try:
    os.mkdir(output)
except Exception as e:
    print(e)
fastafile=open(output+'/'+sample+'.fasta','w')
listfile=open(output+'/'+sample+'.txt','w')
listfile.write('ID\tMUTATION_ID\tSample\tWT.Peptide\tMT.Peptide\tMT.Allele\tWT.Score\tMT.Score\tHLA\tchop_score\n')
num=1
with open(infile,'r')as f:
    for line in f.readlines():
        if 'NeoRank' in line:
            continue
        hla,somaticpep,somaticscore,*_,wildscore,wildpep=line.split('\t')[2:8]
        ch,pos,*_,cds=line.split('\t')[9:13]
        mutation='_'.join([ch,pos,cds,str(len(somaticpep)),hla])
        fastafile.write('>%s|MT|%d|%s\n%s\n'%(sample,num,mutation,somaticpep))
        fastafile.write('>%s|WT|%d|%s\n%s\n'%(sample,num,mutation,wildpep))
        listfile.write('%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n'%(num,mutation,sample,somaticpep,wildpep,hla,wildscore,somaticscore,hla))
        num+=1
listfile.close()
fastafile.close()
