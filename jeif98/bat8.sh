#!/bin/bash
table1-2.sh threaded/intel Angle > table1.pi.tex ; table1-2.sh original/linux Angle > table1.og.tex ; table1-2.sh original/intel Angle > table1.oi.tex ; table1-2.sh threaded/linux Angle > table1.pg.tex
table1-2.sh threaded/intel Value > table2.pi.tex ; table1-2.sh original/linux Value > table2.og.tex ; table1-2.sh original/intel Value > table2.oi.tex ; table1-2.sh threaded/linux Value > table2.pg.tex
