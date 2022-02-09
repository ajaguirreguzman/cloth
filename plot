#!/usr/bin/gnuplot -persist

plot   'energy.out' u 1 w lines tit 'K+V'
replot 'energy.out' u 2 w lines tit 'K'
replot 'energy.out' u 3 w lines tit 'V'
