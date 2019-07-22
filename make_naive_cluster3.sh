#!/bin/bash

FILE="2.3.cluster.txt"

#echo -n "" > $FILE;

for rule in {109..19682}; 
do
    ./build/vouw encode-all -m 2 -b 3 -f 5 -i 2001122102 using -r $rule >> $FILE;
done;
cat $FILE
