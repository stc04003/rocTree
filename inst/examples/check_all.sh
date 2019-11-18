#! /bin/bash

rm .RData
rm .Rhistory
rm *.Rout
rm *.pdf

for f in ex_*
do 
    echo "library(rocTree);source('$f')" > "$f"2.R
    R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-file=$f"out --vanilla < "$f"2.R
    rm "$f"2.R
done
