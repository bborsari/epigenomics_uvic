#!/bin/bash
# get list of non-redundant TSS for each gene contained in the gtf file 
# remember to filter out non protein-coding genes from the gtf in case you don't want them

# define function
# classify-exons 
# (Author: Alessandra Breschi)

classify-exons () {

	cat $1 | awk '
	BEGIN{OFS=FS="\t"; sep=";"}
	$3=="exon" {
		split("", d)
		split($9, a, "; ");
		for(i=1;i<=length(a);i++) {
			split(a[i], b, " ");
			gsub(/"/, "", b[2])
			d[b[1]] = b[2]
		}
		gn_id = d["gene_id"]
		tx_id = d["transcript_id"]
		chr[gn_id] = $1
		strand[gn_id] = $7	
		exon = $1(sep)$4(sep)$5(sep)$7(sep)tx_id
		exon_gn = $1(sep)$4(sep)$5(sep)$7(sep)gn_id


		# Find first and last exons for a transcript
		# ------------------------------------------
		if (min1[tx_id] == "" || min1[tx_id] > $4) {
			min1[tx_id] = $4
			min2[tx_id] = $5
		}

		if (max2[tx_id] == "" || max2[tx_id] < $5) {
			max2[tx_id] = $5
			max1[tx_id] = $4
		}
	
		if (min1[tx_id] == $4) {
			min2[tx_id] = min(min2[tx_id], $5)
		}

		if (max2[tx_id] == $5) {
			max1[tx_id] = max(max1[tx_id], $4)
		}

		ex1[tx_id] = chr[gn_id](sep)min1[tx_id](sep)min2[tx_id](sep)strand[gn_id](sep)tx_id
		ex2[tx_id] = chr[gn_id](sep)max1[tx_id](sep)max2[tx_id](sep)strand[gn_id](sep)tx_id

		# Find constitutive exons
		# -----------------------

		if (txs[gn_id] == "") {
			txs[gn_id] = tx_id
		}

		if (txs[gn_id] !~ tx_id) {
			txs[gn_id] = txs[gn_id]","tx_id
		}

		exonsG[gn_id] = (exonsG[gn_id] == "" ? exon : exonsG[gn_id]","exon)
		exonsT[tx_id] = (exonsT[tx_id] == "" ? exon : exonsT[tx_id]","exon)

		exonCounts[exon_gn] ++
	}

	END {
		for (gn in chr) {

			# Nmber of transcripts for a given gene
			nb_tx = split(txs[gn], t, ",")

			# Iterate over all exons of a gene
		
			for (j=1;j<=nb_tx;j++) {
		
				tx = t[j]
				split(exonsT[tx], e, ",")
	
				# Define first and last exons
				if (strand[gn] == "+") {
					first = ex1[t[j]]
					last = ex2[t[j]]
				}
				if (strand[gn] == "-") {
					first = ex2[t[j]]
					last = ex1[t[j]]
				}

				for(i=1;i<=length(e);i++) {
					pos = "internal"
					if (e[i] == first) {
						pos = "first"
					}
					if (e[i] == last) {
						pos = "last"
					}
					if (length(e) == 1) {
						pos = "unique"
					}
					split(e[i], c, (sep))
					exon_gn = c[1](sep)c[2](sep)c[3](sep)c[4](sep)gn
					const = (exonCounts[exon_gn] == nb_tx ? "constitutive" : "alternative")
					
#					print gn, t[j], e[i], pos, constitutive
					print c[1], c[2], c[3], c[4], pos, const, tx, gn, length(e)
				}
			}
		}
	}

	function min(x,y){return (x<=y) ? x : y;}
	function max(x,y){return (x>=y) ? x : y;}
	'  
}


gtf_file="$1"

classified_exons=$(classify-exons $gtf_file)
for strand in - +; do
	if [ "$strand" == "-" ]; then
		col=3
	else
		col=2
	fi
	grep -F "$strand" <(grep "unique\|first" <(echo "$classified_exons") | awk 'BEGIN{FS=OFS="\t"}{print $1,($2-1),$3,$7,"0",$4,$8}') | sort -u -k7,7 -k"$col","$col"n;
done
