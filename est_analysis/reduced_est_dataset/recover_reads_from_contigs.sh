#! /bin/bash

a=( $( cat contamination_list.csv ) )

for i in {00001..02552}
do
	#echo "Contig$i: "
	for j in "${a[@]}"
	do
    	if [ "contig$i" != "$j" ]
    		then
    			#echo "contig$i != $j"
    			echo "Contig$i: Retained"
				grep contig$i 454ReadStatus.txt | cut -f 1 >> recovered_reads.csv 
			else
				echo "Contig$i: Ignored"
		fi
	done
done

sort recovered_reads.csv -o recovered_reads_sorted.csv

uniq -u recovered_reads_sorted.csv > recovered_reads_uniq.csv