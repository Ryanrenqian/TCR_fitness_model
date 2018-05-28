#!/bin/bash
#=============================================================
# File Name: blastp.sh
# Author: renqian
# mail: renqian@yucebio.com
# Created Time: 2018年05月17日 星期四 10时32分00秒
#=============================================================
#$ -cwd -l vf= -P -q prj.q
export PATH="/home/renqian/bin:/home/renqian/.local/bin:/mnt/cfs/med18b/med/renqian/env/anaconda3/bin:/mnt/cfs/med18b/med/renqian/anaconda3/bin:/mnt/nfs/software/share/jre1.8.0_91/bin:/mnt/nfs/software/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/home/zhoult/perl5/bin:/mnt/cfs/med18b/med/renqian/env/bin:/mnt/nfs/software/opt/texlive/bin/x86_64-linux:/mnt/nfs/user/zhouletian/Or/bin:/mnt/nfs/user/zhouletian/anaconda3/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/opt/dell/srvadmin/bin"
export R_LIBS="/mnt/cfs/med18b/med/renqian/env/lib/R"
export LD_LIBRARY_PATH="/mnt/cfs/med18b/med/renqian/env/lib:/mnt/nfs/software/lib::/opt/MonitorSoftware/lib:/usr/sfw/lib:/usr/local/lib"
sample=$1
output=$2
mkdir -p $output/$sample
/mnt/nfs/software/share/ncbi-blast-2.2.31/bin/blastp \
-query $output/$sample.fasta \
-db /mnt/cfs/med18b/med/renqian/Fitness_model/db/iedb.fasta \
-outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 \
-out $output/$sample/neoantigens_${sample}_iedb.xml
