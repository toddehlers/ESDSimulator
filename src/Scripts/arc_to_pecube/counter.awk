awk '

BEGIN { m=1 }

{

if ( NR == 1 ) { m=$3 }
if ( $3 == m ) { count=count+1 }

}

END { print count }

' temp_sort_file.txt > count.out
