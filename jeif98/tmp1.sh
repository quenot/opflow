#!/bin/bash
paste trc[0-9] | grep -v / > trc
awk '{for(i=1;i<=NF;i++){NUM=NUM?NUM+$i:$i};$(NF+1)=NUM;NUM=""} 1' trc | awk '{printf "%.1f\n", $11/10}'
