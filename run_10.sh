#!/bin/bash
for ((k=1;k<=10;k++)) do
    echo "bin/biddingga > output/apa_sym_2p_5_80_100_1k_$k.csv"
    bin/biddingga > output/apa_sym_2p_5_80_100_1k_$k.csv
done
