#!/bin/bash

rm output

echo "#timing" > timing_864

for((i=1; i<=30; i++))
do

   cp bck_864.box Ar_864.box 

 ./melquiades.x input.dat > output
grep 'Total time in Markov chain :' output | cut -f2- -d":" | cut -c6-21 >> timing_864

echo $i

done
