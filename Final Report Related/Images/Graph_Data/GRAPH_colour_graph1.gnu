reset
set terminal cairolatex pdf color
set output "GRAPH_color_graph4.tex"

set title 'Redshift 10--14'
set xlabel 'f115w-f150w'
set ylabel 'f150w-f200w'
set key top right
#set xrange[0:10]
#set yrange[:0.005]

#f(x) = m*x + c 			#linear
#f(x) = m*exp(-c*x) + d	 	#exponential
#f(x) = m**(c*x) + d 		#power
#f(x) = m*x**(-c) 			#power2
#f(x) = m*log(c*x) + d 		#logarithm

#m = 0.5
#c = 20
#fit f(x) "GRAPH_colour_graphs.txt" using ($1):($3) via m, c, d

#titlef = sprintf("$f(x) = %.3fx + %.3f$", m, c)				#linear
#titlef = sprintf("$f(x) = %.3f\e^{-%.3fx} + %.3f$", m, c, d)		#exponential
#titlef = sprintf("$f(x) = %.3f^{%.3fx} + %.3f$", m, c, d)		#power
#titlef = sprintf("$f(x) = %.3fx^{%.3f}$", m, c)				#power2
#titlef = sprintf("$f(x) = %.3f\log(%.3fx) + %.3f$", m, c, d)	#logarithm

#g(x) = n*log(x) + d
#fit g(x) "StellaDens.txt" using 1:3 via n, d
#titleg = sprintf("%.2flog(x)+%.2f", n, d)

plot\
 "GRAPH_colour_graphs.txt" using 8:7 pt 6 ps 0.7 lw 2 lc rgb "black" notitle

pause -1 "Hit return to continue"

# color plot
#set term post enhanced color solid "Helvetica" 16
# bw plot
#set term post eps enhanced mono dashed "Helvetica" 18

# automatically create pdf file
#set output '| epstopdf --filter --outfile=plot2D.pdf'

# plot again to file:
#replot
