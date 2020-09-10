#!/bin/bash
echo $1
iceg.sh $1 2>&1 >/dev/null | grep real | awk '{print $2}' | sed 's/m/*60+/' | sed 's/s//' | bc -l
