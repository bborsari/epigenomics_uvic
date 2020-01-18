#!/bin/bash
#selects specific columns from a file (columns must be colon-separated)
#the header of the file must have the same length of the rows

[ $# -ge 1 -a -f "$2" ] && input="$2" || input="-"
columns="$1"

awk -v columns="$columns" 'BEGIN{FS=OFS="\t"}{if (NR==1) {n=split ($0, a, "\t");
	k=split (columns, w, ":");
	for (i=1;i<=n;i++){dict_original_file[a[i]]=i};
	for (l=1;l<=k;l++) var=var dict_original_file[w[l]]" "};
	m=split(var, b, " ");
	for (j=1;j<=m;j++)
		{if (j<m) 
			{printf $b[j]"\t"} 
		else if (j==m)
			{printf $b[j]"\n"}}
}' $input
