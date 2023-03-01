set xlabel "x"
set ylabel "y"
m="./data2.txt"
set terminal x11 0
set nokey
set grid
set title 'RFR'
plot m using 1:2 with fsteps