#!/bin/bash
treal6.sh > trc0
treal6.sh > trc1
treal6.sh > trc2
treal6.sh > trc3
treal6.sh > trc4
treal6.sh > trc5
treal6.sh > trc6
treal6.sh > trc7
treal6.sh > trc8
treal6.sh > trc9
paste trc[0-9] | grep -v / | awk '{for(i=1;i<=NF;i++){NUM=NUM?NUM+$i:$i};$(NF+1)=NUM;NUM=""} 1' | awk '{printf "%.1f\n", $11/10}'
