#!/bin/bash

echo Number 1

read numone

echo Number 2

read numtwo

echo Printing from $numone to $numtwo
for (( i=numone; i<numtwo+1; i=i+1))
do
	echo $i
done
