#!/bin/bash

declare -a array=(01 05 1 25)

if [[ ${array[*]} =~ $1 ]]
then 
	echo Cleaning and moving files
	
	./clean_inputs.sh
 	
	cp ../data/Scales/One/* ../data/
 	cp ../data/Fitted/One/Theta_$1/Agg.txt ../data/Fitted/
	cp ../data/Fitted/One/Theta_$1/Theta1.txt ../data/Fitted/
	cp ../data/Fitted/One/Theta_$1/TranParams.csv ../data/	
	exit 0	
else
	echo "Incorrect input!"
	exit 1
fi

