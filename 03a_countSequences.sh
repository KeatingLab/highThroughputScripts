#!/bin/bash
input=$1
output=$2/seq_summary

mkdir $output

for each in `ls $input`;
    do sort $input/$each | uniq -c > $output/$each.smry
done
