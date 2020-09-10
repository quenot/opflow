#!/bin/bash
$1/bin/opflow synthetic/$2/$3 synthetic/$2/$3 cylind -PGM -M 1.5 -C synthetic/$2/vitvue.uwo -CSCALE 0.1 $4 | grep $5 | head -1 | awk '{print $4}'  | sed 's/+/\\pm/' | tr -d '\n'
