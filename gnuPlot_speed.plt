file = "gnuPlot_values.dat"
consts = "gnuPlot_consts.dat"
stats file name "A"
stats file u 8 name "B"
stats consts name "C"

deltaX = C_max_x

set terminal gif animate delay 7 size 1500,1500
set size square
set output "_gnuPlot_speed.gif"

#set palette rgbformula 1,10,1

if (B_min != B_max) {
    set cbrange [B_min:B_max]
}
set xrange [A_min_x:A_max_x]
set yrange [A_min_y:A_max_y]

set cblabel "Length"
set xlabel "Width"
set ylabel "Height"

do for [i=0:int(A_blocks-1)] {
    plot file index i u 1:2:(0.5*$6):(0.5*$7):($8 != 0 ? $8 : NaN) w vectors head size 2*deltaX,15,15 filled lc palette title ""
}