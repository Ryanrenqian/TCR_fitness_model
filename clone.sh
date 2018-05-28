#!/bin/bash
#=============================================================
# File Name: Fitness_model.sh
# Author: renqian
# mail: renqian@yucebio.com
# Created Time: Wed 16 May 2018 02:25:49 PM CST
#=============================================================
#$ -cwd -l vf= -P -q prj.q
export PATH="/home/renqian/bin:/home/renqian/.local/bin:/mnt/cfs/med18b/med/renqian/env/anaconda3/bin:/mnt/nfs/software/share/jre1.8.0_91/bin:/mnt/nfs/software/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/mnt/cfs/med18b/med/renqian/env/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/opt/dell/srvadmin/bin"
export R_LIBS="/mnt/cfs/med18b/med/renqian/env/lib/R"
export LD_LIBRARY_PATH="/mnt/cfs/med18b/med/renqian/env/lib:/mnt/nfs/software/lib::/opt/MonitorSoftware/lib:/usr/sfw/lib:/usr/local/lib"
source activate py27
Normal_BAM=$1
Tumor_BAM=$2
OUTPUT=$3
Sample=$4
VAF=$5
# Step 1 : call SNP
/mnt/cfs/med18b/med/renqian/env/lib/R/facets/extcode/snp-pileup /mnt/nfs/database/hg19/dbsnp/147/All_20160408.vcf.gz $OUTPUT/$Sample/${Sample}.csv.gz $Normal_BAM $Tumor_BAM

# Step 2 : call CNV
Rscript /mnt/cfs/med18b/med/renqian/Fitness_model/facet.R $OUTPUT/$Sample/${Sample} $OUTPUT/$Sample/${Sample}.cnv

# Step3 : prepare data
/mnt/cfs/med18b/med/renqian/Fitness_model/normalize_data.py $OUTPUT/$Sample/${Sample}.cnv.txt $VAF $OUTPUT/$Sample/${Sample}.cnv_input.txt $OUTPUT/$Sample/${sample}.ssm_input.txt

# Step4 : Calculate Clone
/mnt/cfs/med18b/med/renqian/Fitness_model/phylowgs/evolve.py $OUTPUT/$Sample/${Sample}.cnv_input.txt $OUTPUT/$Sample/${sample}.ssm_input.txt

#PS,如果准备从VCF文件中读取突变信息,请在第三步使用phylowgs目录下的paser脚本
#该方法仅支持单样本,如需多样本,请修改第三步脚本,或使用phylowgs目录下的paser脚本(也需要修改)
