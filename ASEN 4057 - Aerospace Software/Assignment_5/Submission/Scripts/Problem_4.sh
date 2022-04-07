#!/bin/bash
# ASEN 4057 HW5 Prob4
# Author: Tyler Candler
# Collabrators: none
# Created: 2/22/2022
# Last Edited: 2/25/2022

#Input: Filename
#Output: Average grade to the console

#cd to the directory with the contained files
cd /home/tyca6175/CANDLER/Assignments/Assignment_5/Submission/Files/txt

#rename input for convenience
file=$1

numbers=$( cut $file -d ';' -f 2 )
#echo $num
#initialize the sum and count
sum=0;
N=0;
#for loops through numbers
for currentnum in $numbers; do
	#sum numbers and add to counter
	let sum=$sum+$currentnum
	let N+=1;
done
#report results to command window
echo "Average for $file is: "
echo "scale=1; $sum / $N" | bc
#echo "Average grade in $file is $avg"

