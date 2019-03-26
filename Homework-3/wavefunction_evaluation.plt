#  File: wavefunction_evaluation.plt 
#
#  GNUplot plot file for wavefunction_evaluation.x output.
#  
#  Programmer:  Daniel Votaw  votaw@nscl.msu.edu
# 
#  Revision history
#   03/26/2019  Taken from wavefunction.plt and modified.

# Put timestamp on plot.
set timestamp

# Titles and labels.
set title 'Wavefunction comparison'
set xlabel 'Dimension N'
set ylabel 'Summed relative deviation'

# Legend placement.
set key right

# Set scales.
set xrange [0:21]
set yrange [0:200]

f(x) = a*x + b
fit [0:21] f(x) "wavefunction_evaluation.dat" using ($1):($2) via a,b 
fit_title = sprintf("%-+4.1f*x %-+4.1f",a,b)

# Set terminal to X11.
set term x11

# Plot exact and approximate wavefunctions.
plot "wavefunction_evaluation.dat" using ($1):($2) title 'Data', \
     a*x + b title fit_title

# Save the image.   
set out "wavefunction_evaluation.ps"
set term postscript color enhanced
replot
