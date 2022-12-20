#!/bin/bash

echo "clean up all test cases for pecube"

if [[ ! $1 ]];
then
    echo "no date/time provided"
    exit 1
fi

echo "clean up for date $1"

for i in case_*;
do
    echo "test case: '$i'"
    cd $i
    rm -v out_$1*.txt
    rm -v core
    cd output
    rm -vrf Pecube_$1*
    cd ..
    cd ..
done

