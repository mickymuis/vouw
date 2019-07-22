#!/bin/bash

FILE="3.2.cluster.txt"

echo -n "" > $FILE;

for rule in {0..255}; 
do
    ./build/vouw encode-all -m 3 -b 2 -f 20 -i 00110 using -r $rule >> $FILE;
done;
cat $FILE
