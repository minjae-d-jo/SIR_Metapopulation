#!/bin/bash

N=$1
L=$2
kappa=$3
p=$4
nRepeat=$5
uid=$6

[[ $# -ne 6 ]] && echo "Error: You need to give 6 arguments." && exit 1


params=${N}_${L}_${kappa}_${p}_${nRepeat}_${uid}
logFile=Log/SIR_$params.table
outputFile=Result/SIR_$params.table


exec 2>> $logFile

echo "$(date +%Y-%m-%d_%H:%M:%S): start $$ $(hostname)" >> $logFile


Bin/SIR_Metapopulation $N $L $kappa $p $nRepeat > $outputFile #2> $outputFile2
ret=$?;
echo -n "$(date +%Y-%m-%d_%H:%M:%S): " >> $logFile
if [ $ret -eq 0 ]; then
	echo "normal exit"
else
	echo "abnormal exit: $ret"
fi >> $logFile

#for((t=1; t<$nRepeat; t++)); do
#	binn/DP_ContactModelSteady $p $T $N #>> $outputFile #2>> $outputFile2
#	#sleep 1
#done >> $outputFile #&2>> $logFile
#echo "$(date +%Y-%m-%d_%H:%M:%S): end $$ $(hostname)" >> $logFile



