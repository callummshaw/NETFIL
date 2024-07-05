#!/bin/bash

declare -a array=(01 05 1 25)

if [[ ${array[*]} =~ $1 ]]
then 
	echo Cleaning and moving files
	
	./clean_inputs.sh
 	
	cp ../data/Scales/Many/* ../data/
 	cp ../data/Fitted/Many/Theta_$1/Agg.txt ../data/Fitted/
	cp ../data/Fitted/Many/Theta_$1/Theta1.txt ../data/Fitted/
	cp ../data/Fitted/Many/Theta_$1/Work.txt ../data/Fitted/
	cp ../data/Fitted/Many/Theta_$1/TranParams.csv ../data/
	exit 0	
else
	echo "Incorrect input!"
	exit 1
fi

