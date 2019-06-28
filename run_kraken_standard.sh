#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
source ~/.bashrc

#soft
kraken='/export4/xielu_group/data175/oy/metagenomics/soft/kraken/kraken'
kraken_mpa_report='/export4/xielu_group/data175/oy/metagenomics/soft/kraken/kraken-mpa-report'
#input
db='/export4/xielu_group/data175/oy/metagenomics/soft/krakendb/standard'

#input
indir='/export4/xielu_group/data175/oy/zhongnan/DNA31/Separate'

#output
outdir='/export4/xielu_group/data175/oy/metagenomics/out'

if [ ! -d $outdir/kraken ]; then
mkdir $outdir/kraken
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

kraken_out_dir=$outdir'/kraken/'$OUT_DIR

echo $kraken_out_dir
if [ ! -d $kraken_out_dir ]; then
mkdir -p $kraken_out_dir
fi


input1=$indir/$FILE1
input2=$indir/$FILE2

kraken_out_dir_out=$kraken_out_dir/$OUT_DIR'.txt'
kraken_out_dir_mpa_report=$kraken_out_dir/$OUT_DIR'.tab.txt'


touch $kraken_out_dir_out
touch $kraken_out_dir_mpa_report


date_start=`date|awk -F"[ :]" '{print $4*3600 + $5*60 +$6}'`

echo "-----------$OUT_DIR kraken start at $date_start ------------------"

echo "$kraken --db $db --threads 40 --fastq-input --gzip-compressed --paired $input1 $input2 > $kraken_out_dir_out"

$kraken --db $db --threads 40 --fastq-input --gzip-compressed --paired $input1 $input2 > $kraken_out_dir_out

echo "$kraken_mpa_report --db $db $kraken_out_dir_out > $kraken_out_dir_mpa_report"

$kraken_mpa_report  --db $db  $kraken_out_dir_out > $kraken_out_dir_mpa_report

date_end=`date|awk -F"[ :]" '{print $4*3600 + $5*60 +$6}'`

echo "-----------$OUT_DIR kraken end at $date_end ------------------"
time=`expr "$date_end" - "$date_start"`
echo "kraken Time takes for $OUT_DIR last for $time seconds"

done < $1

