#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
source ~/.bashrc

/export4/xielu_group/data175/oy/metagenomics/soft/metacv_2_3_0/metacv classify /export4/xielu_group/data175/oy/metagenomics/soft/metacvdb/db/cvk6_2059 /export4/xielu_group/data175/oy/metagenomics/out/bowtie/cmp/cmp.fq.1 /export4/xielu_group/data175/oy/metagenomics/out/bowtie/cmp/cmp.fq.2 /export4/xielu_group/data175/oy/metagenomics/out/metacv/cmp/



