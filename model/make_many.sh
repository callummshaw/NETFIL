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

	if [ $1 == 01 ]
	then
		sed -Ei "2 s/[^,]*/0.1/2" ../data/TranParams.csv
	elif [ $1 == 05 ]
	then
		sed -Ei "2 s/[^,]*/0.5/2" ../data/TranParams.csv   
	elif [ $1 == 25 ]
	then
 		sed -Ei "2 s/[^,]*/2.5/2" ../data/TranParams.csv   
	else 
		sed -Ei "2 s/[^,]*/1/2" ../data/TranParams.csv   
	fi
	exit 0	
else
	echo "Incorrect input!"
	exit 1
fi

