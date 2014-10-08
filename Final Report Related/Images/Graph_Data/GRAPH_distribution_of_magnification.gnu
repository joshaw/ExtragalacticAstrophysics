reset
set terminal cairolatex pdf color
set output "GRAPH_distribution_of_magnification1.tex"

set title ""
set xlabel 'Magnitude ($M$)'
set ylabel 'Number of Galaxies'
set key at 14.8,8
set xrange[9.5:15]
#set yrange[-2.3:-1.2]

#f(x) = m*x + c 			#linear
#f(x) = exp(-m*x + c)	 	#exponential
#f(x) = m**(c*x) + d 		#power
#f(x) = m*x**(-c) 			#power2
#f(x) = m*log(c*x) + d 		#logarithm

#m = 0.5
#c = 20
#fit f(x) "GRAPH_Mag_vs_galaxies_in_time_data.txt" using ($1):($3) via m, c

#titlef = sprintf("$f(x) = %.3fx + %.3f$", m, c)				#linear
#titlef = sprintf("$f(x) = \e^{-%.3fx + %.3f}$", m, c)			#exponential
#titlef = sprintf("$f(x) = %.3f^{%.3fx} + %.3f$", m, c, d)		#power
#titlef = sprintf("$f(x) = %.3fx^{%.3f}$", m, c)				#power2
#titlef = sprintf("$f(x) = %.3f\log(%.3fx) + %.3f$", m, c, d)	#logarithm

#g(x) = n*log(x) + d
#fit g(x) "StellaDens.txt" using 1:3 via n, d
#titleg = sprintf("%.2flog(x)+%.2f", n, d)

plot\
 "GRAPH_distribution_of_magnification_data.txt" using 1:2 pt 6 ps 0.7 lw 4 lc rgb "#4d81be" t '\SI{0.1e6}{\second}', \
 "GRAPH_distribution_of_magnification_data.txt" using 1:3 pt 2 ps 0.7 lw 4 lc rgb "#c05048" t '\SI{0.2e6}{\second}', \
 "GRAPH_distribution_of_magnification_data.txt" using 1:4 pt 1 ps 0.7 lw 4 lc rgb "#9bbb59" t '\SI{0.5e6}{\second}', \
 "GRAPH_distribution_of_magnification_data.txt" using 1:5 pt 1 ps 0.7 lw 4 lc rgb "#71588f" t '\SI{0.5e6}{\second}', \
 "GRAPH_distribution_of_magnification_data.txt" using 1:6 pt 1 ps 0.7 lw 4 lc rgb "#db843d" t '\SI{0.5e6}{\second}'
#with xyerrorbars pt 6 ps 0.7 lw 2 lc rgb "black" notitle, f(x) lw 3 lc rgb "red" t titlef

pause -1 "Hit return to continue"

# color plot
#set term post enhanced color solid "Helvetica" 16
# bw plot
#set term post eps enhanced mono dashed "Helvetica" 18

# automatically create pdf file
#set output '| epstopdf --filter --outfile=plot2D.pdf'

# plot again to file:
#replot
