#!/mnt/cfs/med18b/med/renqian/env/anaconda3/bin/python
###################################################################
# File Name: normalize_data.py
# Author: renqian
# mail: renqian@yucebio.com
# Created Time: Tue 16 Jan 2018 09:25:59 AM CST
#=============================================================
import sys
import numpy as np
print(sys.argv)
cnv,snv,cnr,ssm=sys.argv[1:]
# initial base params
hetsnp_rate=7e-4
avg_ssms_in_tumour = 3000
depths=[]
s=open(ssm,'w')
s.write('id\tgene\ta\td\tmu_r\tmu_v\n')
with open(snv,'r')as f:
    num=0
    gene={}
    mu_r=0.999
    mu_v=0.499
    gene_attr={}
    for line in f.readlines():
        if 'VAF' in line:
            continue
        eles=line.strip().split('\t')
        ref,alt,chm,pos=eles[1:5]
#        print(ref,alt)
        if chm=='X':
            continue
        sid='s'+str(num)
        depths.append(int(alt)+int(ref))
#        if (alt,ref)==('2','2'):
#            continue
        gene[sid]='%s\t%s\t%s\t%d\t%3f\t%3f\n'%(sid,'_'.join((chm,pos)),ref,int(alt)+int(ref),mu_r,mu_v)
        gene_attr.setdefault(chm,{})[pos]=sid
        num+=1
        s.write(gene[sid])
avg_depth=np.nanmedian(depths)
o=open(cnr,'w')
o.write('cnv\ta\td\tssms\tphysical_cnvs\n')
with open(cnv,'r')as f:
    num=0
    for line in f.readlines():
        if  'chromosome' in line or 'NA' in line:
            continue
        chm,st,ed,cn,lcn,mcn,cp=line.strip().split()[1:]
        if (int(mcn),int(lcn)) !=(2,1):
            continue
        if chm in 'XY':
            continue
        # calc a d
        D_max=np.round(avg_ssms_in_tumour * avg_depth).astype(np.int)
        D = int(np.round((int(ed) - int(st)) * hetsnp_rate*avg_depth))
        d=np.minimum(D_max,D)
        a=int((1 - float(cp)/2) * d)
        ssms=[]
        chmdict=gene_attr.get(chm,{})
        for key in chmdict.keys():
            if int(st)<= int(pos)<=int(ed):
                ssms.append(','.join(([gene_attr[chm][key],lcn,mcn])))
#                s.write(gene[gene_attr[chm][key]])
        cid='c'+str(num)
        num+=1
        phycnv='chrom=%s,start=%s,end=%s,major_cn=%s,minor_cn=%s,cell_prev=%3f'%(chm,st,ed,mcn,lcn,float(cp))
        o.write('%s\t%d\t%d\t%s\t%s\n'%(cid,a,d,';'.join(ssms),phycnv))
s.close()
o.close()

