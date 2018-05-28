#!/bin/bash
#=============================================================
# File Name: neoantigen_fitness.sh
# Author: renqian
# mail: renqian@yucebio.com
# Created Time: 2018年05月17日 星期四 10时01分31秒
#=============================================================
#$ -cwd -l vf= -P -q prj.q
export PATH="/home/renqian/bin:/home/renqian/.local/bin:/mnt/cfs/med18b/med/renqian/env/anaconda3/bin:/mnt/cfs/med18b/med/renqian/anaconda3/bin:/mnt/nfs/software/share/jre1.8.0_91/bin:/mnt/nfs/software/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/home/zhoult/perl5/bin:/mnt/cfs/med18b/med/renqian/env/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/mnt/nfs/user/zhouletian/Or/bin:/mnt/nfs/user/zhouletian/anaconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/opt/dell/srvadmin/bin"
export R_LIBS="/mnt/cfs/med18b/med/renqian/env/lib/R"
export LD_LIBRARY_PATH="/mnt/cfs/med18b/med/renqian/env/lib:/mnt/nfs/software/lib::/opt/MonitorSoftware/lib:/usr/sfw/lib:/usr/local/lib"
source activate py27
set -x
SAMPLE=$1
OUTPUT=$2
NEOANTIGEN=$3
# STEP1 : 生成多肽
mkdir -p $OUTPUT/$SAMPLE
/mnt/cfs/med18b/med/renqian/Fitness_model/neoantigen_report2data.py $NEOANTIGEN $OUTPUT $SAMPLE && \
# STEP2 ：比对多肽
/mnt/nfs/software/share/ncbi-blast-2.2.31/bin/blastp \
-query $OUTPUT/$SAMPLE.fasta \
-db /mnt/cfs/med18b/med/renqian/Fitness_model/db/iedb.fasta \
-outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 \
-out $OUTPUT/$SAMPLE/neoantigens_${SAMPLE}_iedb.xml
# STEP3 ：计算Fitness
a=26.
k=4.86936
python /mnt/cfs/med18b/med/renqian/Fitness_model/fitness/src/main.py \
$OUTPUT/$SAMPLE.txt \
$OUTPUT/$SAMPLE \
$a $k \
$OUTPUT/neoantigen_fitness_${SAMPLE}.txt
