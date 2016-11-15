#!/bin/bash
#for ((k=300; k<=; k=k+1 ))
#for ((k=27; k<=32; k=k+1 ))  
#do
#echo $k
#./makeAna2.csh $k
#done
rm -rf *dat
./makeAna3.csh 800
rm -rf *dat
./makeAna3.csh 1000
rm -rf *dat
./makeAna3.csh 1400
rm -rf *dat
./makeAna3.csh 2000
rm -rf *dat
./makeAna3.csh 3000
rm -rf *dat
./makeAna3.csh merge
