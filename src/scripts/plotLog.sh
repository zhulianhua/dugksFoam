#!/bin/bash
grep -e 'Temperature'  log |cut -d '=' -f 2 |nl >tlog.dat
grep -e '^Density'  log |cut -d '=' -f 2 |nl >rlog.dat
grep -e '^Velocity'  log |cut -d '=' -f 2 |nl >ulog.dat
gnuplot <<- EOF
    set xlabel "t"
    set ylabel "log(Error)"
    set logscale y
    #set format y "%s*10^{%S}"
    set format y "%.2e"
    set term png
    set output "log.png"
    plot "tlog.dat" using 1:2 with l title "rho Convergence rate", \
         "rlog.dat" using 1:2 with l title "T   Convergence rate", \
         "ulog.dat" using 1:2 with l title "U   Convergence rate"
EOF
rm *log.dat
