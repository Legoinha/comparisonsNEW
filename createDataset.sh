#!/bin/bash

for year in 8
do 
    for num in `seq 1 99` 
    do
        root -q -l 'createDataset.cc('${year}',4,0,100,'${num}')'
        #root -q -l 'createDataset.cc('${year}',4,1,4,'${num}')'
    done
        #root -q -l 'createDataset.cc('${year}',6,1,1,0)'
        #root -q -l 'createDataset.cc('${year}',6,0,1,0)'
done
