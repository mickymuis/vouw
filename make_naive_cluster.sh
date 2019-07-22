#!/bin/bash

FILE="2.2.cluster.txt"

echo -n "" > $FILE;

for rule in {0..15}; 
do
    ./build/vouw encode-all -m 2 -b 2 -f 30 -i 00110 using -r $rule >> $FILE;
done;
cat $FILE
