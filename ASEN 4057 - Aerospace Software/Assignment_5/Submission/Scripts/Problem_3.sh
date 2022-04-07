#!/bin/bash
# ASEN 4057 HW5 Prob3
# Author: Tyler Candler
# Collabrators: none
# Created: 2/22/2022
# Last Edited: 2/25/2022

#Input: a base file name and extension
#Output: Renames all files to a numbered list with the baseneame and file extension provided

name=$1
ext=$2
#cd into directory where pictures are located
cd /home/tyca6175/CANDLER/Assignments/Assignment_5/Submission/Files/jpg

#starting number for naming
count=1
#for loop goes through files with the extension specified by the inputes
for currentfile in *;
do
	#print which file is being renamed	
	echo "$currentfile renamed to $name$count.$ext"

	#rename the current file in the for loop
	mv "$currentfile" "$name$count.$ext"
	#add 1 to the count
	let count+=1;
done


