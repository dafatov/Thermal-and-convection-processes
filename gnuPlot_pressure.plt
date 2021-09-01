file = "gnuPlot_values.dat"
stats file name "A"
stats file u 5 name "B"

set terminal gif animate delay 7 size 1500,1500
set size square
set output "_gnuPlot_pressure.gif"

#set palette rgbformula 33,13,10

if (B_min != B_max) {
    set cbrange [B_min:B_max]
}
set xrange [A_min_x:A_max_x]
set yrange [A_min_y:A_max_y]

set cblabel "Pressure"
set xlabel "Width"
set ylabel "Height"

do for [i=0:int(A_blocks-1)] { plot file index i u 1:2:5 w image}