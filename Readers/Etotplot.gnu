set encoding iso_8859_15
#set terminal postscript enhanced solid color "Helvetica" 20
#set output "gnuplot.ps"
set key 
set yrange [-545.2:-545.35]
set label "Equation of state"
set xlabel "a.u."
set ylabel "E"
set label "Ecal" 
plot   "outfile" u 1:2 w l title "Ecal" lw 2

