#!/bin/bash
#
animals=("a cat" "a dog" "a fish")
for i in "${animals[@]}"
do
	echo $i
done
echo ${animals[0]}
