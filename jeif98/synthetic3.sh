#!/bin/bash
echo -n '$'
bash synthetic1.sh $1 $2 $3 "-WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 6" $4
echo -n '$ & $'
bash synthetic1.sh $1 $2 $3 "-WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 4 -P 2 -N 2 -MF 4 -SUM" $4
echo -n '$ & $'
bash synthetic1.sh $1 $2 $3 "-WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 4 -P 2 -N 2 -MF 4" $4
echo '$ \\'
