#!/bin/bash
time $1/bin/opflow real/iceg real/iceg icegc00 -PGM -WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 6 -M 41.5
time $1/bin/opflow real/iceg real/iceg icegc00 -PGM -WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 4 -MF 4 -M 41.5 -P 2 -N 2
time $1/bin/opflow real/iceg real/iceg icegc00 -PGM -WIN 0.7 -SPW 0.7 -GOA -SMO 4.0 -EPS 3.0 -MB 4 -MF 4 -SUM -M 41.5 -P 2 -N 2
