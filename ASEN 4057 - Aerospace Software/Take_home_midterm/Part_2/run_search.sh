#!bin/bash

#for test

#rsync -a $1 ~/CANDLER/Assignments/Take_home_midterm/Part_2/
cd $1

for currentfile in *.in; do #loop though all files with .in extension
   #echo $currentfile
   #search for requested numbers on each file
   ~/CANDLER/Assignments/Take_home_midterm/Part_2/Program2dot1 $currentfile 17
   ~/CANDLER/Assignments/Take_home_midterm/Part_2/Program2dot1 $currentfile 56
	#./Program2dot1 $currentfile 17
   ~/CANDLER/Assignments/Take_home_midterm/Part_2/Program2dot1 $currentfile 58
  done
