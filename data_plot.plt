set xlabel "x"
set ylabel "y"
m="./data.txt"
set nokey
set grid
set title 'Trees'
plot m using 1:2 with points