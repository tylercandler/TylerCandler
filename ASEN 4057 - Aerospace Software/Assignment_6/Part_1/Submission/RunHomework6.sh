#!bin/bash
#0, 10, 100, 1000, 5000, 10000, 50000, 100000


echo "Testing for objective $1 moon-spacecraft clearance of $2 meters:"

./Homework6 $1 $2 $3


tolerance=$(echo "scale=0; $3 * 10" |bc)
x=$(printf "%0.f" $tolerance)

#echo "$x"

mv ThreeBody_data.csv Optimum_$1_$2_0p$x.csv
