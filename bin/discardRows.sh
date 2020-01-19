#!/bin/bash
#selects the rows from file2 that are not present in file1
#files must not have headers

#usage is:
#discardRows.sh file1 file2

file1="$1"
file2="$2"
awk 'NR==FNR{a[$0];next}!($1 in a)' $file1 $file2
