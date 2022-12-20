#!/bin/bash
# Takes input file and creates formatted text file to be read into Pecube as
# topography file

clear

# Asks for data file and output filename
echo 'Enter ascii three column elevation filename'
read file
echo 'Enter output text file name without extension'
read output

# Creates temporary file
temp='temp_sort_file.txt'

# Checks for input  and output files
if [ -e $file ]
	then echo 'Processing..'
	else
	     echo 'File '$file' does not exist'
	     exit 1
fi

if [ -e $output.txt ]
	then echo 'Output file already exists, choose different name'
	     exit 1
	else echo ''
fi

# Copies contents of the input file to temporary file
# Deletes the first line (header line)
cat $file > $temp
vi -c :1,1d -c :wq $temp

# Converts temporary file to Unix format
dos2unix $temp 

# Runs script counter.awk which finds how many elevations there are for each y value
./counter.awk

# Sets count to the ouput of counter.awk, counts how many lines are in the whole
# file, divides these two values to find how many iterations needed to sort the data
count=`head -1 count.out`
numline=`wc -l $temp | awk '{ print $1 }'`
numobject=`expr $numline / $count`

# while loop that takes each set of values for each y value
# Starts from bottom and works up, it flips the order and then deletes the x and y
# values to just give elevations in the final output file
i=`expr 1`
while [ $i -le $numobject ]
do
	current=$(($count*$i))
	`tail -$current $temp | head -$count | awk '{ print $1 }' >> $output.txt`
	i=`expr $i + 1`
done


# Calulates the spacing between x and y coordinates
long1=`awk 'BEGIN {NR=0}; {if(NR == 1) print $2}' $temp`
long2=`awk 'BEGIN {NR=0}; {if(NR == 2) print $2}' $temp`

next=`expr $count + 1`
lat1=`awk 'BEGIN {NR=0}; {if(NR == '$count') print $3}' $temp`
lat2=`awk 'BEGIN {NR=0}; {if(NR == '$next') print $3}' $temp`

long_space=$(echo "scale=4; $long2-$long1" | bc)
lat_space=$(echo "scale=4; $lat1-$lat2" | bc)

# Prints points and spacing for longitude and latitude
echo 'Number of longitude points: '$count
echo 'Number of latitude points: '$numobject
echo ''
echo 'Spacing of longitude points: '$long_space
echo 'Spacing of latitude points: '$lat_space


rm $temp
rm count.out
	
exit 0
