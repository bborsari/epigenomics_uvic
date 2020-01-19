#!/bin/bash


#********
# USAGE *
#********

display_usage() { 
	echo -e "DESCRIPTION: provided a list of ChIP-seq experiments, returns a matrix of values of the ChIP-seq target for each gene in each experiment \n"
	echo -e "\t--bw <path_to_bw_files> (i.e. a .txt file with the path to the bw files to be used by bwtool, sorted by the desired order of experiments)\n"
	echo -e "\t--bedfile <bedFile> (i.e. gene bodies, TSSs or enhancers)\n"
        echo -e	"\t--target <ChIP-seq target>\n"
	echo -e "\t--outFolder <output folder> (the folder where to save the bwtool summary files; default: cwd)\n"
	echo -e "\t--outFile <output file> (default: target.matrix.tsv)\n"
	echo -e "\t--signal <max/mean/median/sum> (i.e. the column to retrieve from the bwtool summary file; default: median)\n"
	echo -e "\t--peaks <path_to_peakCalling_output> (i.e. a .txt file with the path to the peak calling bed files, sorted by the same order of bw file; not required for marks of type 'broad')\n"
	echo -e "\t--keep <yes/no> (whether to keep or not the temporary files; default: no)\n"
	echo -e "\tIMPORTANT: bw and peaks .txt files must have the same length and the same order of experiments.\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi


if [  $# -le 3  ]
then
	echo -e "ERROR: insufficient number of arguments\n"
    	display_usage
        exit 1
fi



#******************
# READING OPTIONS *
#******************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in

        --bedfile)
        bedFile="$2"
        shift
        ;;
    	
	--bw)
    	bw="$2"
    	shift
    	;;
    	
	--outFolder)
    	outFolder="$2"
    	shift
    	;;

	--outFile)
	outFile="$2"
	shift
	;;

	--target)
	target="$2"
	shift
	;;

	--peaks)
	peaks="$2"
	shift
	;;

	--signal)
	signal="$2"
	shift
	;;

	--keep)
	keep="$2"
	;;
	*)
	
	;;
	esac
	shift
done


: ${outFolder:="."}
: ${outFile:="$target"".matrix.tsv"}
: ${keep:="no"}
: ${signal:="median"}


if [[ "$signal" != "max" && "$signal" != "mean" && "$signal" != "median" && "$signal" != "sum" ]]
then
	echo "ERROR: signal must be of max, mean, median or sum"
	display_usage
	exit 1
fi


echo "Reading options .."
echo "path to bw files =" "${bw}"
echo "bed file =" "${bedFile}"
echo "output folder =" "${outFolder}"
echo "output file =" "${outFile}"
echo "targeted histone mark =" "${target}"
echo "path to peak calling output =" "${peaks}"
echo "signal =" "${signal}"
echo -e "keep tmp files =" "${keep}\n"


#*******************
# DEFINE FUNCTIONS *
#*******************

# 1. define function
# add_header

add_header () {

	[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

	awk 'BEGIN{FS=OFS="\t"}{if(NR==1){
					for(i=1;i<=(NF-2);i++)
					{printf "V"i"\t"};
					printf "V"NF-1"\n";
					print $0}
				else
					{print $0}}' $input
}


# 2. define function
# selectColumns

selectColumns () {

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

}


#********
# BEGIN *
#********

#********
# 1. Perform intersection between the bedfile with the regions of interest
# and the bedfile with peaks for each experiment

# In case a region doesn't intersect with any peak
# it returns the region itself, otherwise it returns
# the coordinates of the intersection between the region and the peak
#********

echo -e "retrieve coordinates of overlapping regions..\n"

n=0
cat $peaks | while read file; do
	n=$(echo $n+1 | bc -l)
		
	bedtools intersect -a $bedFile -b $file -wao | \
	awk 'BEGIN { FS=OFS="\t" }{if ( $NF==0 ) {print $1, $2, $3, $7"_"$4, $5, $6 } else {start = ($2>$9)?$2:$9; end = ($3<$10)?$3:$10; print $1, start, end, $7"_"$4, $5, $6 }}' > "$outFolder"/"$(basename $file .bed)".coordinates.intersection.bed
		
		
	if [[ "$n" == 1 ]]
	then		
		echo "$outFolder"/$(basename $file .bed)".coordinates.intersection.bed" > "$outFolder"/"$target".path.coordinates.txt
	else
		echo "$outFolder"/$(basename $file .bed)".coordinates.intersection.bed" >> "$outFolder"/"$target".path.coordinates.txt
	fi
done

	
if [[ $(awk 'END{print NR}' $peaks) != $(awk 'END{print NR}' "$outFolder"/"$target".path.coordinates.txt) ]]
then
	echo "ERROR: length of coordinates file is different from length of input peak calling file"
	exit 1
fi


#********
# 2. Compute bwtool summary stats
# (mean, median, max, sum signal)
# over the set of selected regions 
# for each experiment
#********

echo -e "bwtool summary of overlapping regions..\n"

n=$(awk 'END{print NR}' "$outFolder"/"$target".path.coordinates.txt)
for i in `seq 1 $n`; do
	bed_coordinates=$(head -$i "$outFolder"/"$target".path.coordinates.txt | tail -1)
	file=$(head -$i $bw | tail -1)
	bwtool summary -header -with-sum -keep-bed $bed_coordinates $file "$outFolder"/"$(basename $file .bigWig)".all.genes.summary.tsv
	if [[ "$i" == 1 ]]
	then 
		echo "$outFolder/"$(basename $file .bigWig)".all.genes.summary.tsv" > "$outFolder"/"$target".path.bwtool.summary.txt
	else
		echo "$outFolder/"$(basename $file .bigWig)".all.genes.summary.tsv" >> "$outFolder"/"$target".path.bwtool.summary.txt
	fi
done


#********
# 3. Generate matrix of target's signal
# across genes (rows) and experiments
# (columns). If multiple regions are present
# for a given gene, it will select
# the one with highest signal
#********

echo -e "generate matrix..\n"

file1=$(head -1 "$outFolder"/"$target".path.bwtool.summary.txt)
selectColumns name:"$signal" "$file1" | tail -n+2 |awk 'BEGIN{FS=OFS="\t"}{split($1, a, "_"); print a[1], $2}' | sort -k1,1 -k2,2gr | sort -u -k1,1 > "$outFile"


if [[ "$keep" == "no" ]]
then
	rm $file1
fi

tail -n+2 "$outFolder"/"$target".path.bwtool.summary.txt | while read file; do
	
	selectColumns name:"$signal" "$file" | tail -n+2 | awk 'BEGIN{FS=OFS="\t"}{split($1, a, "_"); print a[1], $2}' | sort -k1,1 -k2,2gr | sort -u -k1,1 > tmp
	if [[ $(diff <(cut -f1 "$outFile") <(cut -f1 tmp)) == "" ]]
	then
		paste "$outFile" <(cut -f2 tmp) > tmp2
		mv tmp2 "$outFile"
		rm tmp
		if [[ "$keep" == "no" ]]
		then
			rm $file
		fi
	else
		echo "ERROR: matrix and tmp files are not sorted in the same way"
		exit 1
	fi

done



#********
# 4. Replace NAs with zeros
#********

echo -e "remove NAs..\n"

sed 's/NA/0/g' "$outFile" > tmp
mv tmp "$outFile"


#*********
# 5. Add header
#*********
	
# echo -e "headerize matrix..\n"

add_header "$outFile" > tmp
mv tmp "$outFile"


#*********
# 6. Remove tmp files
# if desired	
#*********

if [[ "$keep" == "no" ]]
then
	echo -e "remove temporary files..\n"
	cat "$target".path.coordinates.txt | while read file; do
		if [ -e $file ]
		then
			rm $file
		fi
	done
	rm "$target".path.coordinates.txt
fi

if [[ "$keep" == "no" ]]
then
	rm "$target".path.bwtool.summary.txt
fi

