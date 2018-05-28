import re,sys,logging
from optparse import OptionParser
def getopt():
    parser=OptionParser()
    parser.add_option('-t','--tree',dest='tree',help='PhyloWGS tree file')
    parser.add_option('-s','--ssm',dest='ssm',help='SSMs file for building PhyloWGS tree')
    parser.add_option('-o','--output',dest='output',help='outputfile')
    parser.add_option('-v','--verbose',action='store_true',dest='verbose',help='show log or not')
    parser.add_option('-l','--lohhla',dest='lohhla',default=None,help='HLA LOH prediction file')
    parser.add_option('-f','--fitness',dest='fitness',help='fitness of A, R and AXR of each neoantigen')
    parser.add_option('-c','--cutoff',dest='cutoff',type='float',default=0.01,help='cutoff value of HLA LOH Pvalue')
    return parser.parse_args()

def _parse_trees(file):
    with open(file,'r') as f:
        num=-1
        trees={}
        for i in f.readlines():
            if 'id' in i:
                num+=1
                tree={}
                trees[num]=tree
            ele=i.strip().split(',')
            if ele[0].isdigit() and int(ele[3])>0:
                tree[ele[0]]=ele[4].strip().split('; ')
    return trees

def _check_tree(tree):
    totallength=0
    totaleles=[]
    for branch in tree.values():
        totallength+=len(branch)
        totaleles+=branch
    return totallength==len(set(totaleles))

import pandas as pd
import re

if __name__=='__main__':
    options,arg=getopt()
    output=options.output
    logging.basicConfig(filename=output+'.log',level=logging.DEBUG, format=' %(asctime)s - %(levelname)s - %(message)s')
    kp=[]
    allhla=[]
    if not options.verbose:
        logging.disable(logging.DEBUG)
    if  options.lohhla:
        with open(options.lohhla,'r')as f:
            for line in f.readlines():
                if 'region' in line:
                    continue
                ele=line.strip().split()
                lost,keep=ele[-4:-2]
                allhla+=[lost,keep]
                pvalue=ele[-8]
                kp.append(keep)
                if pvalue == 'NA' or float(pvalue)>options.cutoff:
                    kp.append(lost)
        logging.debug('keep:\t%s'%str(set(kp)))
        logging.debug('lost:\t%s'%str(set(allhla)-set(kp)))
        kp=['_'.join(i.split('_')[0:4]) for i in kp]
#        ls=[lost.split('_')[0:4] for lost in set(allhla)-set(kp)]
        ls=set(['_'.join(lost.split('_')[0:4]) for lost in set(allhla)-set(kp)])-set(kp)
    else:
        ls=None
#尽管生成多棵树，只取第一个结果输出
    trees=_parse_trees(options.tree)
    for tree in trees.values():
        if _check_tree(tree):
            break
    logging.debug('number_of_branch:\t%d'%(len(tree)))
# 注释sid
    sidinfo={}
    with open(options.ssm,'r')as f:
        for i in f.readlines():
            if 'id' in i:
                continue
            sid,mut=i.split()[:2]
            sidinfo[mut]=sid
    #使用pandas操作表
    fitness=pd.read_table(options.fitness)
    fitness['HLA']=fitness.apply(lambda x:re.sub('\W','_',x[1].split('_')[-1].lower()),axis=1)
    def getClone(x):
        maxfit={}
        mut='_'.join(x[1].split('_')[:2])
        if sidinfo.get(mut,None):
            for branch in tree.keys():
                if sidinfo[mut] in tree[branch]:
                    return branch
        return None
    fitness['Clonal']=fitness.apply(getClone,axis=1)
    maxfit={}
    for i in fitness.index:
        clonal=fitness.loc[i]['Clonal']
        if clonal==None:
            continue
        maxfit[clonal]=max(fitness.loc[i]['NeoantigenRecognitionPotential'],maxfit.get(clonal,0))
    fitness['Clonal_Fitness']=fitness.apply(lambda x: maxfit.get(x['Clonal'],x['NeoantigenRecognitionPotential']),axis=1)
    if ls:
    # 0 represent keep, 1 represent loss
        fitness['LOHHLA']=fitness.apply(lambda x: 1 if x["HLA"] in ls else 0,axis=1)
        logging.debug('neoantigen of LOHHLA:\t%d'%sum(fitness['LOHHLA']))
        logging.debug('clonal  of LOHHLA:\t%d'%sum(fitness[fitness['LOHHLA']==1]['Clonal']=='1'))
    print(-sum(maxfit.values()))
    fitness.to_csv(output+'.txt',index=False,sep='\t')
    with open(output+'sum.txt','w')as o:
        o.write('%s\n'%str(-sum(maxfit.values())))

# 计算样本的fitness
#    with open(options.fitness,'r')as f:
#        maxfit={}
#        for i in f.readlines():
#            mut='_'.join(i.split()[1].split('_')[:2])
#            hla=i.split()[1].split('_')[-1]
#            hla=re.sub('\W','_',hla.lower())
#            logging.debug('HLA:%s'%hla)
#            if kp:
#                if hla not in kp:
#                    continue
#            logging.debug('sid:%s'%sidinfo.get(mut,''))
#            if sidinfo.get(mut,None):
#                for branch in tree.keys():
#                    if  sidinfo[mut] in tree[branch]:
#                        logging.debug('fit:%s'%i.strip().split()[-1])
#                        clonal=branch
#                        maxfit[branch]=max(maxfit.get(branch,0),float(i.strip().split()[-1]))
#        logging.debug('mutation: %s'%mut)
#    fit=-sum(maxfit.values())
#    o=open(output,'w')
#    o.write('Sample Fitness: %f\n'%fit)
#    o.close()
#    print(output,fit)
#    logging.debug('sample fitness:%s'%str(maxfit))
#    logging.debug('branch %s'%str(tree[branch]))
