#!/bin/bash

echo "comparing two ages files for all results"

if [[ -z $1 ]];
then
    echo "you must provide a date and time, for ex.:"
    echo "$0 2013_10_28__11_04 2013_12_05__09_28"
    exit 1
fi

if [[ -z $2 ]];
then
    echo "you must provide a date and time, for ex.:"
    echo "$0 2013_10_28__11_04 2013_12_05__09_28"
    exit 1
fi

for i in case_*;
do
    echo "test case: '$i'"

    FILE1=""
    FILE2=""

    if [[ -f "$i/output/Pecube_$1/Ages_tec001.dat" ]];
    then
        FILE1="$i/output/Pecube_$1/Ages_tec001.dat"
    elif [[ -f "$i/output/Pecube_$1/Ages_tec002.dat" ]];
    then
        FILE1="$i/output/Pecube_$1/Ages_tec002.dat"
    fi

    if [[ -f "$i/output/Pecube_$2/Ages_tec001.dat" ]];
    then
        FILE2="$i/output/Pecube_$2/Ages_tec001.dat"
    elif [[ -f "$i/output/Pecube_$2/Ages_tec002.dat" ]];
    then
        FILE2="$i/output/Pecube_$2/Ages_tec002.dat"
    fi

    echo "comparing these two files:"
    echo $FILE1
    echo $FILE2

    diff $FILE1 $FILE2

done

