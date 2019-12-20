#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
source ~/.bashrc

#soft
metacv="/export4/xielu_group/data175/oy/metagenomics/soft/metacv_2_3_0/metacv"

#db
metacvdb='/export4/xielu_group/data175/oy/metagenomics/soft/metacvdb/db/cvk6_2059'

#input
indir='/export4/xielu_group/data175/xxjian/DNA_22_20190805/samples'

#output
outdir='/export4/xielu_group/data175/xxjian/DNA_22_20190805'


if [ ! -d $outdir/metacv ]; then
mkdir $outdir/metacv
fi


while read  LINE
do
echo $LINE
echo $OUT_DIR
FILE1=`echo $LINE | cut -d \; -f 1`
FILE2=`echo $LINE | cut -d \; -f 2`
OUT_DIR=`echo $FILE1 | cut -d \/ -f 1`
echo $FILE1
echo 'OUT_DIR='$OUT_DIR

metacv_out_dir=$outdir'/metacv/'$OUT_DIR

echo $metacv_out_dir
if [ ! -d $metacv_out_dir ]; then
mkdir -p $metacv_out_dir
fi

input1=$indir/$FILE1
input2=$indir/$FILE2

#gzip -d $input1
#gzip -d $input2


date_start=`date|awk -F"[ :]" '{print $4*3600 + $5*60 +$6}'`
echo "-----------$OUT_DIR kraken start at $date_start -----------------"

echo "$metacv classify $metacvdb $input1 $input2 $metacv_out_dir --threads=64"
# $metacv classify /export4/xielu_group/data175/oy/metagenomics/soft/metacvdb/db/cvk6_2059 $input1 $input2 $metacv_out_dir
$metacv classify $metacvdb $input1 $input2 $metacv_out_dir --threads=64

date_end=`date|awk -F"[ :]" '{print $4*3600 + $5*60 +$6}'`
echo "-----------$OUT_DIR metacv end at $date_end ------------------"
time=`expr "$date_end" - "$date_start"`
echo "metacv Time takes for $OUT_DIR last for $time seconds"

done <$1
