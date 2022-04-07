#!/bin/bash
# ASEN 4057 HW5 Prob2
# Author: Tyler Candler
# Collabrators: none
# Created: 2/22/2022
# Last Edited: 2/22/2022

#Input: A Name
#Output: Whether the input is a file or directory and the permissions associated with it


echo What is the name of the file/directory?

read varname
#determine if the input is a file
if [ -e $varname ];  then
#file exists
echo "$varname exists"
	#check if directory
	if [ -d $varname ]; then
		echo "$varname is a directory"
	fi
	#check if regular file
	if [ -f $varname ]; then
		echo "$varname is a regular file"
	fi
#set permissions to a variable
permissions=$(ls -l $varname| head -c 10 | cut -c 2-4)
#check permissions variable for wrx and output the proper permissions on the command line
echo %permissions
	if [ $permissions = *r* ]; then
		echo "User can read";
	fi
	if [ $permissions = *w* ]; then
                echo "User can read";
        fi
	if [ $permissions = *x* ]; then
                echo "User can read";
  	fi
	else
	echo "Error: file does not exist or is mispelled"
fi



