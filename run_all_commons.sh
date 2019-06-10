#!/bin/bash
for ((i=1;i<=5;i++)); 
do
    for ((j=i;j<=5;j++)); 
    do
	for ((k=j;k<=5;k++)); 
	do
	    for ((m=k;m<=5;m++)); 
	    do
		echo "bin/common_value $i $j $k $m > output/common_second_$i$j$k$m.csv"
		bin/common_value $i $j $k $m > output/common_second_$i$j$k$m.csv
	    done
	done
    done
done
