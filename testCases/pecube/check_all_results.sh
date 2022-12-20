#!/bin/bash

# general tests for all cases:

function check_for_output_file {
    # check if output file exists
    if [[ ! -f $1 ]];
    then
        echo "FAILED: '$2', logfile '$1' does not exist"
        exit 1
    fi
}

function check_if_exit_success {
    # check if pecube finished normally or if it crashed
    # if it finished normally, then the run time information is given at the end of the file

    NUM_OF_LINES=$(tail -n 3 $1 | grep " Total times :" | wc -l)

    if [[ $NUM_OF_LINES != 3 ]];
    then
        echo "FAILED: '$2', logfile '$1' does not contain time information at the end"
        exit 1
    fi
}

function check_for_nan {
    # check for NaN in any output file

    grep -I -i "NaN" $1/*
    EXIT_STATUS=$?

    if [[ $EXIT_STATUS -eq 0 ]];
    then
        echo "FAILED: '$2', found NaN in folder '$1'"
        exit 1
    fi
}

function check_for_correct_ages {
    # check that all ages are in correct order
    #
    # Tecplot header:
    # 1 AHe Age - Farley, 2000 (Ma)
    # 2 AHe Age - a= 20.0 um (Ma)
    # 3 AHe Age - a= 40.0 um (Ma)
    # 4 AHe Age - a= 70.0 um (Ma)
    # 5 AHe Age - eU= 12.4 ppm (Ma)
    # 6 AHe Age - eU= 61.8 ppm (Ma)
    # 7 AHe Age - eU=247.0 ppm (Ma)
    # 8 AFT Age (Ma)
    # 9 KfeldAr Age (Ma)
    # 10 ZHe Age (Ma)
    # 11 ZFT Age (Ma)
    # 12 MuscAr Age (Ma)
    # 13 BioAr Age (Ma)
    # 14 BioAr2 (Ma)
    # 15 MAr Age (Ma)
    # 16 Ap U-Th / Pb (Ma)
    # 17 HornAr Age (Ma)


    TEMP_FILE=$(mktemp)

    for age_file in $1/Ages_tec*.dat;
    do
        tail -n+5 $age_file | awk 'NF == 20 {$1=$2=$3=""; print}' > $TEMP_FILE
        awk '$1 > $8 {printf("line %d: AHe %f > AFT %f\n", NR + 4, $1, $8); exit 1}' $TEMP_FILE
        EXIT_STATUS=$?

        awk '$8 - $9 > 2.0 {printf("line %d: AFT %f > KfeldAr %f\n", NR + 4, $8, $9); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$9 > $10 {printf("line %d: KfeldAr %f > ZHe %f\n", NR + 4, $9, $10); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$10 > $11 {printf("line %d: ZHe %f > ZFT %f\n", NR + 4, $10, $11); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$11 - $12 > 1.0 {printf("line %d: ZFT %f > MuscAr %f\n", NR + 4, $11, $12); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$12 > $13 {printf("line %d: MuscAr %f > BioAr %f\n", NR + 4, $12, $13); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$13 > $15 {printf("line %d: BioAr %f > MAr %f\n", NR + 4, $13, $15); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$13 > $16 {printf("line %d: BioAr %f > Ap U-Th / Pb %f\n", NR + 4, $13, $16); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$13 > $17 {printf("line %d: BioAr %f > HornAr %f\n", NR + 4, $13, $17); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        awk '$15 > $16 {printf("line %d: MAr %f > Ap U-Th / Pb %f\n", NR + 4, $15, $16); exit 1}' $TEMP_FILE
        (( EXIT_STATUS |= $? ))

        if [[ $EXIT_STATUS -eq 1 ]];
        then
            echo "FAILED: '$2', ages not in correct order in file '$age_file'"
            exit 1
        fi
    done

    rm $TEMP_FILE

}

echo "checking all results for pecube"

if [[ -z $1 ]];
then
    echo "you must provide a date and time, for ex. 2013_10_28__11_04"
    exit 1
fi


for i in case_*;
do
    echo "test case: '$i'"
    cd $i

    OUTPUT_FILE=out_$1.txt
    OUTPUT_DIR=output/Pecube_$1

    # run general tests:
    check_for_output_file $OUTPUT_FILE $i
    check_if_exit_success $OUTPUT_FILE $i
    check_for_nan $OUTPUT_DIR $i
    check_for_correct_ages $OUTPUT_DIR $i

    # run specific tests:
    ./check_results.sh $1 $i

    cd ..
done

echo "all tests passed without an error"

