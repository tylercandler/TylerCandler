#!/bin/bash
# ASEN 4057 HW5 Prob4
# Author: Tyler Candler
# Collabrators: none
# Created: 2/22/2022
# Last Edited: 2/25/2022

#Input: Compressed Tar File
#Output: Extracted and sorted directory
tarfile=$1

#cd to related directoyr
cd /home/tyca6175/CANDLER/Assignments/Assignment_5/Submission/Files

#create new directory for new files
mkdir EXTRACTED_CLEANED

#unzip files and move to newly created directory
tar -xf $tarfile -C EXTRACTED_CLEANED
#cd to new directory
cd EXTRACTED_CLEANED
cd MISC

#for loop through all files in directory
for file in *; do
	#check for file extension
	ext="${file##*.}"
	#does file have extension
	if [[ $ext == txt ]]; then
		if [ -d $ext ];then
			mv $file $ext
		fi
		else 
			mkdir $ext
			mv $file $ext
	
	fi
	
	
	
done

#cd back to previous directory
cd /home/tyca6175/CANDLER/Assignments/Assignment_5/Submission/Files
#compress directory
tar -cvf MISC_COMPRESSED_CLEAN.tar.gz EXTRACTED_CLEANED

