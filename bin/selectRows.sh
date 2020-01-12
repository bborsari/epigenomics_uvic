#!/bin/bash
#selects the rows from file1 that are present in file2
#if the files have different headers, it does not include the header in the output (to be solved)

file1="$1"
file2="$2"
awk 'NR==FNR{a[$0];next} $1 in a' $file1 $file2

