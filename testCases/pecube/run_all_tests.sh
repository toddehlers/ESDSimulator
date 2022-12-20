#!/bin/bash

echo "running all test cases for pecube"

for i in case_*;
do
    echo "test case: '$i'"
    cd $i
    ./run_test.sh &
    cd ..
done

